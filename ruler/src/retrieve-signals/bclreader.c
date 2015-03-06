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


struct BCLData *
load_bcl_file(const char *filename)
{
    struct BCLData *bcl=NULL;
    uint32_t nclusters, i;
    uint8_t *base, *quality;
    FILE *fp;

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "load_bcl_file: Can't open file %s.\n", filename);
        return NULL;
    }

    if (fread(&nclusters, sizeof(nclusters), 1, fp) < 1) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        goto onError;
    }

    nclusters = le32toh(nclusters);
    bcl = malloc(sizeof(struct BCLData));
    if (bcl == NULL) {
        perror("load_bcl_file");
        goto onError;
    }

    bcl->nclusters = nclusters;
    base = bcl->base = malloc(nclusters);
    quality = bcl->quality = malloc(nclusters);
    if (base == NULL || quality == NULL) {
        perror("load_bcl_file");
        goto onError;
    }

    if (fread(base, nclusters, 1, fp) < 1) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        printf("nclusters=%u file=%s\n", nclusters, filename);
        goto onError;
    }

    fclose(fp);

    for (i = 0; i < nclusters; i++, base++, quality++) {
        if (*base == 0) {
            *quality = NOCALL_QUALITY + PHRED_BASE;
            *base = NOCALL_BASE;
        }
        else {
            *quality = (*base >> 2) + PHRED_BASE;
            *base = CALL_BASES[*base & 3];
        }
    }

    return bcl;

  onError:
    if (bcl != NULL) {
        if (bcl->base != NULL)
            free(bcl->base);
        if (bcl->quality != NULL)
            free(bcl->quality);
        free(bcl);
    }

    fclose(fp);

    return NULL;
}


void
free_bcl_data(struct BCLData *data)
{
    free(data->base);
    free(data->quality);
    free(data);
}


void
format_basecalls(char *seq, char *qual, struct BCLData **basecalls,
                 int ncycles, uint32_t clusterno)
{
    uint32_t i;

    for (i = 0; i < ncycles; i++) {
        *seq++ = basecalls[i]->base[clusterno];
        *qual++ = basecalls[i]->quality[clusterno];
    }

    *seq = *qual = 0;
}
