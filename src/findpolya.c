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
#include "tailseq-sigproc.h"


uint32_t
find_polya(const char *seq, size_t seqlen, struct PolyAFinderParameters *params)
{
    int max_term_mod, *weights;
    int longest_i, longest_j, longest_length;
    int i, j;

    max_term_mod = params->max_terminal_modifications;
    weights = params->weights;

    if (seqlen < max_term_mod)
        max_term_mod = seqlen;

    longest_i = longest_j = seqlen + 1;
    longest_length = -1;

    /* calculate match scores for all possible [i, j] */
    for (i = 0; i < max_term_mod; i++) {
        int curlength, scoresum;

        scoresum = weights[(int)seq[i]];
        if (longest_length < 1 && scoresum > 0) {
            longest_length = 1;
            longest_i = longest_j = i;
        }

        for (j = i + 1, curlength = 2; j < seqlen; j++, curlength++) {
            scoresum += weights[(int)seq[j]];

            if (scoresum > 0 && curlength > longest_length) {
                longest_i = i;
                longest_j = j;
                longest_length = curlength;
            }
        }
    }

    if (longest_length < 0)
        return 0;

    while (seq[longest_i] != 'T' && longest_i <= longest_j)
        longest_i++;

    while (seq[longest_j] != 'T' && longest_j >= longest_i)
        longest_j--;

    if (longest_j + 1 - longest_i < params->min_polya_length)
        return 0;
    else
        return ((uint32_t)longest_i << 16) | (longest_j + 1);
}
