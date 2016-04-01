/*
 * signalproc.c
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
#include <limits.h>
#include <endian.h>
#include "tailseq-sigproc.h"


int
measure_polya_len(struct IntensitySet *intensities, const char *seq,
                  int ncycles, struct PolyARulerParameters *params)
{
    return -1;
}


int
load_color_matrix(float *mtx, const char *filename)
{
    float origmtx[NUM_CHANNELS * NUM_CHANNELS];
    FILE *fp;
    int i;
    
    fp = fopen(filename, "rt");
    if (fp == NULL) {
        fprintf(stderr, "Color matrix <%s> could not be opened.\n", filename);
        return -1;
    }

    for (i = 0; i < NUM_CHANNELS * NUM_CHANNELS; i++) {
        if (fscanf(fp, "%f", origmtx + i) < 1) {
            fclose(fp);

            fprintf(stderr, "Color matrix <%s> could not be loaded.\n", filename);
            return -1;
        }
    }

    fclose(fp);

    if (inverse_4x4_matrix(origmtx, mtx) < 0) {
        fprintf(stderr, "Color matrix could not be inverted.\n");
        return -1;
    }

    return 0;    
}

