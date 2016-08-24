/*
 * findpolya.c
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
#include "tailseq-import.h"


uint32_t
find_polya(const char *seq, size_t seqlen, struct PolyAFinderParameters *params)
{
    short *polyA_weights, *nonA_weights;
    int max_term_mod, nonA_score;
    int best_i, best_j, best_length, best_score;
    int i, j, min_polya_len;

    max_term_mod = params->max_terminal_modifications;
    polyA_weights = params->weights_polyA;
    nonA_weights = params->weights_nonA;
    min_polya_len = params->min_polya_length;

    if (seqlen < max_term_mod)
        max_term_mod = seqlen;

    best_i = best_j = -1;
    best_length = best_score = -1;
    nonA_score = 0;

    /* calculate match scores for all possible [i, j] */
    for (i = 0; i < max_term_mod; i++) {
        int curlength, scoresum;

        scoresum = nonA_score + polyA_weights[(int)seq[i]];
        if (min_polya_len <= 1 && scoresum > best_score) {
            best_i = best_j = i;
            best_length = 1;
            best_score = scoresum;
        }

        for (j = i + 1, curlength = 2; j < seqlen; j++, curlength++) {
            scoresum += polyA_weights[(int)seq[j]];

            if (scoresum > best_score && curlength >= min_polya_len) {
                best_i = i;
                best_j = j;
                best_length = curlength;
                best_score = scoresum;
            }
        }

        nonA_score += nonA_weights[(int)seq[i]];
    }

    if (best_length < (int)params->min_polya_length || best_score < 1)
        return 0;
    else
        return ((uint32_t)best_i << 16) | (best_j + 1);
}
