/*
 * bclreader.c
 *
 * Copyright (c) 2015 Hyeshik Chang
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
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "tailseq-retrieve-signals.h"


#define NOCALL_QUALITY      2
#define NOCALL_BASE         'N'
#define PHRED_BASE          33
#define BATCH_BLOCK_SIZE    65536

static const char CALL_BASES[4] = "ACGT";


struct BCLReader *
open_bcl_file(const char *filename)
{
    struct BCLReader *bcl;
    FILE *fp;

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "load_bcl_file: Can't open file %s.\n", filename);
        return NULL;
    }

    bcl = malloc(sizeof(struct BCLReader));
    if (bcl == NULL) {
        fclose(fp);
        return NULL;
    }

    if (fread(&bcl->nclusters, sizeof(bcl->nclusters), 1, fp) < 1) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        fclose(fp);
        free(bcl);
        return NULL;
    }

    bcl->nclusters = le32toh(bcl->nclusters);
    bcl->fptr = fp;
    bcl->read = 0;

    return bcl;
}


void
close_bcl_file(struct BCLReader *bcl)
{
    if (bcl == BCLREADER_OVERRIDDEN)
        return;

    if (bcl->fptr == NULL)
        fclose(bcl->fptr);
    free(bcl);
}


int
load_bcl_data(struct BCLReader *bcl, struct BCLData *data, uint32_t nclusters)
{
    uint32_t toread;

    if (bcl->read >= bcl->nclusters) {
        data->nclusters = 0;
        return 0;
    }

    if (bcl->read + nclusters >= bcl->nclusters) /* does file has enough clusters to read? */
        toread = bcl->nclusters - bcl->read; /* all clusters left */
    else
        toread = nclusters;

    if (fread(data->basequality, toread, 1, bcl->fptr) < 1) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        return -1;
    }

    data->nclusters = toread;
    bcl->read += toread;

    return 0;
}


struct BCLData *
new_bcl_data(uint32_t size)
{
    struct BCLData *bcl;

    bcl = malloc(sizeof(size_t) + size);
    if (bcl == NULL) {
        perror("new_bcl_data");
        return NULL;
    }

    bcl->nclusters = 0;

    return bcl;
}


void
free_bcl_data(struct BCLData *data)
{
    free(data);
}


void
format_basecalls(char *seq, char *qual, struct BCLData **basecalls,
                 int ncycles, uint32_t clusterno)
{
    uint32_t i;

    for (i = 0; i < ncycles; i++) {
        uint8_t bq=basecalls[i]->basequality[clusterno];

        if (bq == 0) {
            *qual++ = NOCALL_QUALITY + PHRED_BASE;
            *seq++ = NOCALL_BASE;
        }
        else {
            *qual++ = (bq >> 2) + PHRED_BASE;
            *seq++ = CALL_BASES[bq & 3];
        }
    }

    *seq = *qual = 0;
}


struct BCLReader **
open_bcl_readers(const char *msgprefix, const char *datadir, int lane, int tile, int ncycles,
                 struct AlternativeCallInfo *altcalls)
{
    int16_t cycleno, bcllaststart;
    char path[PATH_MAX];
    struct BCLReader **readers;

    readers = malloc(sizeof(struct BCLReader *) * ncycles);
    if (readers == NULL) {
        perror("open_bcl_readers");
        return NULL;
    }

    memset(readers, 0, sizeof(struct BCLReader *) * ncycles);

    /* mark overridden cycles by alternative calls not to open a reader for them. */
    for (; altcalls != NULL; altcalls = altcalls->next) {
        int16_t last_cycle=altcalls->first_cycle + altcalls->reader->ncycles;

        for (cycleno = altcalls->first_cycle; cycleno < last_cycle; cycleno++)
            readers[cycleno] = BCLREADER_OVERRIDDEN;
    }

    bcllaststart = -1;

    for (cycleno = 0; cycleno < ncycles; cycleno++) {
        if (readers[cycleno] != NULL) {
            if (bcllaststart >= 0) {
                printf("%sUsing BCL base calls for cycle %d-%d.\n", msgprefix,
                        bcllaststart+1, cycleno);
                bcllaststart = -1; 
            }

            continue;
        }

        snprintf(path, PATH_MAX, "%s/BaseCalls/L%03d/C%d.1/s_%d_%04d.bcl", datadir, lane,
                 cycleno+1, lane, tile);

        readers[cycleno] = open_bcl_file(path);
        if (readers[cycleno] == NULL) {
            while (--cycleno >= 0)
                close_bcl_file(readers[cycleno]);
            return NULL;
        }

        if (bcllaststart < 0)
            bcllaststart = cycleno;
    }

    if (bcllaststart >= 0)
        printf("%sUsing BCL base calls for cycle %d-%d.\n", msgprefix,
                bcllaststart+1, cycleno);

    return readers;
}

void
close_bcl_readers(struct BCLReader **readers, int ncycles)
{
    int cycleno;

    for (cycleno = 0; cycleno < ncycles; cycleno++)
        close_bcl_file(readers[cycleno]);
}
