/*
 * polyaruler.c
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
#include <string.h>
#include <stdint.h>
#include <zlib.h>
#include <math.h>
#include <assert.h>
#include "../sigproc-flags.h"


#define SCORE_LINEBUF_SIZE      16384
#define MAX_NUM_CYCLES          1024
#define TAGINFO_LINEBUF_SIZE    2048


struct SignalRecordHeader {
    uint32_t clusterno;
    int16_t first_cycle;            /* left-most cycle number of polyA start */
    int16_t valid_cycle_count;      /* cycle count with valid signals */
};

static inline int
min_int(int a, int b)
{
    if (a < b) return a;
    else return b;
}

static float *
load_score_cutoffs(const char *filename, const char *tileid,
                   ssize_t *ncycles)
{
    FILE *fp;
    char linebuf[SCORE_LINEBUF_SIZE], *bptr, *cutofftok;
    size_t tileid_len;
    float *cutoffs;
    int i;

    tileid_len = strlen(tileid);
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s.\n", filename);
        return NULL;
    }

    /* Search for the line for the designated tile. */
    for (;;) {
        if (fgets(linebuf, sizeof(linebuf), fp) == NULL) {
            fprintf(stderr, "No poly(A) score cutoffs found for tile %s.\n",
                    tileid);
            fclose(fp);
            return NULL;
        }

        if (strlen(linebuf) >= tileid_len &&
                memcmp(linebuf, tileid, tileid_len) == 0 &&
                linebuf[tileid_len] == '\t')
            break;
    }

    fclose(fp);

    cutoffs = malloc(sizeof(float) * MAX_NUM_CYCLES);
    if (cutoffs == NULL)
        return NULL;

    bptr = linebuf + tileid_len + 1;
    for (i = 0; (cutofftok = strsep(&bptr, "\t\r\n, ")) != NULL; i++) {
        if (*cutofftok == 0)
            break;

        if (i >= MAX_NUM_CYCLES) {
            fprintf(stderr, "Exceeded the maximum allowed number of cycles. "
                            "Adjust MAX_NUM_CYCLES in " __FILE__ ".\n");
            free(cutoffs);
            return NULL;
        }

        cutoffs[i] = atof(cutofftok);
    }

    *ncycles = i;

    return cutoffs;
}

static int
measure_polya_length(const float *read_scores, ssize_t read_num_cycles,
                     int first_cycle,
                     const float *score_cutoffs, ssize_t cutoffs_num_cycles)
{
    ssize_t i, physical_cycle;
    int score_cumsum, score_cumsum_max, score_argmax_cycle;

    score_cumsum = score_cumsum_max = 0;
    score_argmax_cycle = -1;
    physical_cycle = first_cycle;

    for (i = 0; i < read_num_cycles && physical_cycle < cutoffs_num_cycles;
            i++, physical_cycle++) {
        float cutoff = score_cutoffs[physical_cycle];

        if (isnan(cutoff))
            continue;

        score_cumsum += -1 + (read_scores[i] >= cutoff) * 2;
        if (score_cumsum > score_cumsum_max) {
            score_cumsum_max = score_cumsum;
            score_argmax_cycle = i;
        }
    }

    return score_argmax_cycle;
}

static int16_t *
process_polya_ruling(const char *filename, const float *score_cutoffs,
                     ssize_t cutoffs_num_cycles, int minimum_polya_len,
                     size_t *ret_elements)
{
    gzFile fp;
    ssize_t record_size, bytesread;
    int16_t *polya_measurements;
    struct {
        uint32_t elemsize;
        uint32_t total_clusters;
        uint32_t max_cycles;
    } header;

    fp = gzopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open the input file: %s\n", filename);
        return NULL;
    }

    if (gzread(fp, &header, sizeof(header)) != sizeof(header)) {
        fprintf(stderr, "Failed to read the header from %s.\n", filename);
        gzclose(fp);
        return NULL;
    }

    if (header.elemsize != sizeof(float)) {
        fprintf(stderr, "The file %s was written in a machine with "
                        "different architecture.\n", filename);
        gzclose(fp);
        return NULL;
    }

    polya_measurements = malloc(sizeof(int16_t) * header.total_clusters);
    if (polya_measurements == NULL) {
        perror("process_polya_ruling");
        gzclose(fp);
        return NULL;
    }

    memset(polya_measurements, 0xff,
           sizeof(int16_t) * header.total_clusters); /* fill with -1 */
    record_size = header.max_cycles * sizeof(float);

    for (;;) {
        int polya_len;
        struct TSRecord {
            struct SignalRecordHeader header;
            float scores[header.max_cycles];
        } rec;

        bytesread = gzread(fp, (void *)&rec, sizeof(struct TSRecord));
        if (bytesread == 0)
            break;

        if (bytesread < record_size) {
            fprintf(stderr, "Unexpected end of file.\n");
            return NULL;
        }

        polya_len = measure_polya_length(rec.scores,
                rec.header.valid_cycle_count, rec.header.first_cycle,
                score_cutoffs, cutoffs_num_cycles);

        if (polya_len >= minimum_polya_len)
            polya_measurements[rec.header.clusterno] = polya_len;
    }

    gzclose(fp);
    *ret_elements = header.total_clusters;

    return polya_measurements;
}

static int
output_corrected_polya_measurements(const char *taginfo_file,
                                    int16_t *polya_measurements,
                                    size_t nclusters,
                                    const char *tile_id)
{
    gzFile fp;

    fp = gzopen(taginfo_file, "rt");
    if (fp == NULL)
        return -1;

    for (;;) {
        char linebuf[TAGINFO_LINEBUF_SIZE];
        char *bptr, *token, *mods, *umi;
        uint32_t clusterno, flags;
        int i;

        if (gzgets(fp, linebuf, TAGINFO_LINEBUF_SIZE) == NULL) {
            int errno;
            (void)gzerror(fp, &errno);
            if (errno == 0)
                break; /* end-of-file */
            else {
                fprintf(stderr, "Error occurred on reading %s.\n",
                        taginfo_file);
                gzclose(fp);
                return -1;
            }
        }

        bptr = linebuf;
        token = strsep(&bptr, "\t\n\r");
        if (token == NULL)
            continue;

        clusterno = atoi(token);
        if (polya_measurements[clusterno] < 0) {
            /* Poly(A) length is not revised. Bypass the line. */
            token[strlen(token)] = '\t';
            fputs(tile_id, stdout);
            fputc('\t', stdout);
            fputs(linebuf, stdout);
            continue;
        }

        mods = umi = NULL;
        flags = 0;

        for (i = 0; (token = strsep(&bptr, "\t\n\r")) != NULL; i++)
            switch (i) {
            case 0: flags = atoi(token); break;
            case 1: /*polya_prelim = atoi(token);*/ break;
            case 2: mods = token; break;
            case 3: umi = token; break;
            default: break;
            }

        flags |= PAFLAG_MEASURED_FROM_FLUORESCENCE;
        printf("%s\t%d\t%d\t%d\t%s\t%s\n", tile_id, clusterno, flags,
               polya_measurements[clusterno], mods, umi);
    }

    return 0;
}

int
main(int argc, char *argv[])
{
    ssize_t cutoffs_num_cycles;
    float *score_cutoffs;
    int16_t *polya_measurements;
    size_t nclusters;
    const char *tile_id, *signals_file, *cutoffs_file, *taginfo_file;
    int minimum_polya_len;

    if (argc < 6) {
        fprintf(stderr, "Usage: %s {tile id} {signals} {cutoffs} {min polya} "
                        "{taginfo}\n", argv[0]);
        return 1;
    }

    tile_id = argv[1];
    signals_file = argv[2];
    cutoffs_file = argv[3];
    minimum_polya_len = atoi(argv[4]);
    taginfo_file = argv[5];

    /* Load per-cycle poly(A) score cutoffs table */
    score_cutoffs = load_score_cutoffs(cutoffs_file, tile_id, &cutoffs_num_cycles);
    if (score_cutoffs == NULL)
        return 1;

    /* Measure the lengths of poly(A) tails */
    polya_measurements = process_polya_ruling(signals_file, score_cutoffs,
            cutoffs_num_cycles, minimum_polya_len, &nclusters);
    free(score_cutoffs);
    if (polya_measurements == NULL)
        return 2;

    /* Apply the measurements to the existing taginfo */
    if (output_corrected_polya_measurements(taginfo_file,
                polya_measurements, nclusters, tile_id) < 0) {
        free(polya_measurements);
        return 3;
    }

    free(polya_measurements);

    return 0;
}
