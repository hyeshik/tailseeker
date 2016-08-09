/*
 * tailseq-dedup-perfect.c
 *
 * Copyright (c) 2016 Hyeshik Chang
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * - Hyeshik Chang <hyeshik@snu.ac.kr>
 */

#define _BSD_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "../sigproc-flags.h"

#define DEFAULT_BUFFER_SIZE     1024
#define MAX_TILENAME_LEN        63
#define MAX_MODIFICATION_LEN    50
#define MAX_UMI_LEN             50
#define MAX_LINE_LEN            511
#define BUFFER_EXPANSION_FACTOR 1.5

struct TagInfo {
    char line[MAX_LINE_LEN+1];

    char tilename[MAX_TILENAME_LEN+1];
    int clusterno;
    int flags;
    int polyA_len;

    char modifications[MAX_MODIFICATION_LEN+1];
    size_t modification_len;

    char umi[MAX_UMI_LEN+1];
    size_t umi_len;
};

struct TagInfoQueue { /* a circular queue */
    struct TagInfo *el;

    ssize_t size;
    ssize_t front;
    ssize_t rear;

    int num_duplicates;
    int highest_priority;

    FILE *statsout;
};


#define tqueue_length(q)    (((q)->front + (q)->size - (q)->rear) % (q)->size)
#define tqueue_head(q)      ((q)->el[(q)->front])
#define tqueue_tail(q)      ((q)->el[(q)->rear])
#define tqueue_next(q, eln) ((eln + 1) % (q)->size)
#define tqueue_inc(q)       do { (q)->front = tqueue_next((q), (q)->front); } while(0)
#define tqueue_dec(q)       do { (q)->rear = tqueue_next((q), (q)->rear); } while(0)
#define tqueue_isempty(q)   ((q)->front == (q)->rear)
#define tqueue_isfull(q)    (tqueue_length(q) == (q)->size - 1)
#define tqueue_foreach(q, v) for ((v) = (q)->rear; (v) != (q)->front; (v) = tqueue_next((q), (v)))
#define tqueue_empty(q)     do { (q)->rear = (q)->front; } while(0)


static inline int
calculate_tag_prority(int flags)
{
#define SET(value, flagbit) (((value) & (flagbit)) != 0)
#define NOTSET(value, flagbit) (((value) & (flagbit)) == 0)
    return
        NOTSET(flags, PAFLAG_BARCODE_HAS_MISMATCHES) * 1024 +
        NOTSET(flags, PAFLAG_DELIMITER_NOT_FOUND) * 512 +
        NOTSET(flags, PAFLAG_DELIMITER_HAS_MISMATCH) * 256 +
        SET(flags, PAFLAG_MEASURED_FROM_FLUORESCENCE) * 128 +
        NOTSET(flags, PAFLAG_DARKCYCLE_OVER_THRESHOLD) * 64 +
        NOTSET(flags, PAFLAG_BALANCER_BIASED) * 32 +
        NOTSET(flags, PAFLAG_BALANCER_SIGNAL_BAD) * 16 +
        NOTSET(flags, PAFLAG_DARKCYCLE_EXISTS) * 8 +
        NOTSET(flags, PAFLAG_DELIMITER_IS_SHIFTED) * 4 +
        NOTSET(flags, PAFLAG_NO_POLYA_DETECTED) * 2 +
        NOTSET(flags, PAFLAG_MEASURED_USING_NAIVE_RULER);
#undef SET
#undef NOTSET
}

static struct TagInfoQueue *
tqueue_new(ssize_t size)
{
    struct TagInfoQueue *queue;

    queue = malloc(sizeof(struct TagInfoQueue));
    if (queue == NULL)
        return NULL;

    queue->el = calloc(size, sizeof(struct TagInfo));
    if (queue->el == NULL) {
        free(queue);
        return NULL;
    }

    queue->size = size;
    queue->front = 0;
    queue->rear = 0;
    queue->num_duplicates = 0;
    queue->highest_priority = -1;
    queue->statsout = NULL;

    return queue;
}

static void
tqueue_free(struct TagInfoQueue *queue)
{
    free(queue->el);
    if (queue->statsout != NULL)
        fclose(queue->statsout);
    free(queue);
}

static int
tqueue_expand(struct TagInfoQueue *queue, float factor)
{
    ssize_t newsize, wptr, rptr;
    struct TagInfo *newel;

    newsize = (ssize_t)(queue->size * factor);
    newel = calloc(newsize, sizeof(struct TagInfo));
    if (newel == NULL)
        return -1;

    wptr = 0;
    tqueue_foreach(queue, rptr)
        memcpy(&newel[wptr++], &queue->el[rptr], sizeof(struct TagInfo));

    assert(wptr == tqueue_length(queue));

    free(queue->el);
    queue->el = newel;
    queue->size = newsize;
    queue->rear = 0;
    queue->front = wptr;

    return 0;
}

static int
tqueue_append(struct TagInfoQueue *queue)
{
    /* The new entry to be appended is passed as located at the front pointer. */
    struct TagInfo *current;
    int tagprio;

    current = &tqueue_head(queue);
    tagprio = calculate_tag_prority(current->flags);

    if (tqueue_isempty(queue) || tagprio > queue->highest_priority) {
        tqueue_empty(queue);
        tqueue_inc(queue);
        queue->highest_priority = tagprio;
        queue->num_duplicates++;
        return 0;
    }

    if (tqueue_isfull(queue) &&
            tqueue_expand(queue, BUFFER_EXPANSION_FACTOR) < 0)
        return -1;

    queue->num_duplicates++;
    if (tagprio >= queue->highest_priority) {
        tqueue_inc(queue);
        queue->highest_priority = tagprio;
    }

    return 0;
}

static inline int
parse_line(struct TagInfo *taginfo)
{
    size_t length;
    char *pos, *start;

    /* Locate tilename and copy */
    pos = strchr(taginfo->line, '\t');
    if (pos == NULL)
        return -1;
    length = (size_t)(pos - taginfo->line);
    memcpy(taginfo->tilename, taginfo->line, length);
    taginfo->tilename[length] = '\0';

    /* Locate cluster number and parse it */
    start = pos + 1;
    taginfo->clusterno = strtoul(start, &pos, 10);
    if (pos == NULL)
        return -1;

    /* Locate flags and parse it */
    start = pos + 1;
    taginfo->flags = strtol(start, &pos, 10);
    if (pos == NULL)
        return -1;

    /* Locate poly(A) length and parse it */
    start = pos + 1;
    taginfo->polyA_len = strtol(start, &pos, 10);
    if (pos == NULL)
        return -1;

    /* Locate 3' end modification sequence and copy */
    start = pos + 1;
    pos = strchr(start, '\t');
    if (pos == NULL)
        return -1;
    length = (size_t)(pos - start);
    if (length > 0)
        memcpy(taginfo->modifications, start, length);
    taginfo->modifications[length] = '\0';
    taginfo->modification_len = length;

    /* Trim the UMI from the input buffer as we'll append a new column to the lines. */
    *pos = '\0';

    /* Locate UMI sequence and copy */
    start = pos + 1;
    pos = strchr(start, '\n');
    if (pos == NULL)
        return -1;
    length = (size_t)(pos - start);
    if (length > 0)
        memcpy(taginfo->umi, start, length);
    taginfo->umi[length] = '\0';
    taginfo->umi_len = length;

    return 0;
}


static int
process_tag_duplicates(struct TagInfoQueue *queue)
{
    ssize_t ptr, length;

    length = tqueue_length(queue);
    if (length == 0)
        return 0;

    if (length == 1) {
        printf("%s\t%d\n", tqueue_tail(queue).line, queue->num_duplicates);
        if (queue->statsout != NULL) {
            struct TagInfo *taginfo=&tqueue_tail(queue);
            fprintf(queue->statsout, "%s\t%d\t%d\t%d\n", taginfo->tilename, taginfo->clusterno,
                    taginfo->polyA_len, taginfo->polyA_len);
        }
    }
    else {
        int polyA_len_sum, nearest_ptr;
        float mean_polyA_len, nearest_dist;

        /* Calculate the arithmetic mean of poly(A) lengths of duplicates. */
        polyA_len_sum = 0;
        tqueue_foreach(queue, ptr) {
            polyA_len_sum += queue->el[ptr].polyA_len;
        }
        mean_polyA_len = (float)polyA_len_sum / length;

        /* Find the tag instance that is nearest to the mean. */
        nearest_ptr = queue->rear;
        nearest_dist = fabsf(queue->el[queue->rear].polyA_len - mean_polyA_len);
        for (ptr = tqueue_next(queue, queue->rear);
                ptr != queue->front; ptr = tqueue_next(queue, ptr)) {
            float dist=fabsf(queue->el[ptr].polyA_len - mean_polyA_len);
            if (dist < nearest_dist) {
                nearest_ptr = ptr;
                nearest_dist = dist;
            }
        }

        /* Output the representative tag with the mean poly(A) length. */
        {
            struct TagInfo *rep;
            int final_polyA;

            rep = &queue->el[nearest_ptr];
            final_polyA = rep->polyA_len >= 0 ? (int)roundf(mean_polyA_len) : -1;
            printf("%s\t%d\t%d\t%d\t%s\t%d\n",
                rep->tilename, rep->clusterno, rep->flags, final_polyA,
                rep->modifications, queue->num_duplicates);

            /* Output poly(A) length calls for accuracy assessments of clones */
            if (queue->statsout != NULL) {
                fprintf(queue->statsout, "%s\t%d\t%d", rep->tilename, rep->clusterno,
                        final_polyA);

                for (ptr = tqueue_next(queue, queue->rear);
                        ptr != queue->front; ptr = tqueue_next(queue, ptr))
                    fprintf(queue->statsout, "\t%d", queue->el[ptr].polyA_len);
                fprintf(queue->statsout, "\n");
            }
        }
    }

    return 0;
}


int
main(int argc, char *argv[])
{
    struct TagInfoQueue *queue;
    int lineno;

    queue = tqueue_new(DEFAULT_BUFFER_SIZE);
    if (queue == NULL) {
        perror("tqueue_new");
        goto onError;
    }

    if (argc >= 2) {
        queue->statsout = fopen(argv[1], "w");
        if (queue->statsout == NULL) {
            fprintf(stderr, "ERROR: Cannot open %s.", argv[1]);
            goto onError;
        }
    }

    for (lineno = 0;
         fgets(tqueue_head(queue).line, MAX_LINE_LEN, stdin) != NULL;
         lineno++) {
        struct TagInfo *current;

        current = &tqueue_head(queue);

        if (parse_line(current) < 0) {
            fprintf(stderr, "Could line parse line %d: %s", lineno, current->line);
            break;
        }

        if (!tqueue_isempty(queue) &&
                strcmp(current->umi, tqueue_tail(queue).umi) != 0) {
            if (process_tag_duplicates(queue) < 0) {
                perror("process_tag_duplicates");
                goto onError;
            }

            queue->num_duplicates = 0;
            queue->highest_priority = -1;
        }

        if (tqueue_append(queue) < 0) {
            perror("tqueue_append");
            goto onError;
        }
    }

    if (process_tag_duplicates(queue) < 0) {
        perror("process_tag_duplicates");
        goto onError;
    }

    tqueue_free(queue);

    if (ferror(stdin)) {
        perror("fgets");
        return -1;
    }

    return 0;

  onError:
    if (queue != NULL)
        tqueue_free(queue);

    return -1;
}
