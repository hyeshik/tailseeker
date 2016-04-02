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
#include <assert.h>
#include <math.h>
#include "tailseq-sigproc.h"


static void
decrosstalk_intensity(float *out, const struct IntensitySet *original,
                      const float *colormatrix)
{
    int j;

    for (j = 0; j < NUM_CHANNELS; j++) {
#if NUM_CHANNELS == 4
        out[j] =
            colormatrix[0] * original->value[0] +
            colormatrix[1] * original->value[1] +
            colormatrix[2] * original->value[2] +
            colormatrix[3] * original->value[3];
#else
  #error The channel count is not supported yet
#endif
        colormatrix += NUM_CHANNELS;
    }
}


static int
check_balance_minimum(const char *seq, const struct PolyARulerParameters *params)
{
    short found[256];
    char *base;
    int i;

    memset(found, 0, sizeof(short) * 256);

    for (i = params->balancer_start; i < params->balancer_end; i++) {
        printf("%c", seq[i]);
        found[(int)seq[i]]++;
    }

    for (base = "ACGT"; *base != '\0'; base++)
        if (found[(int)*base] < params->balancer_minimum_occurrence)
            return -1;

    return 0;
}


static int
probe_signal_ranges(float *signal_range_high, float *signal_range_low,
                    const struct IntensitySet *intensities,
                    const struct PolyARulerParameters *params)
{
    int npos=params->balancer_num_positive_samples;
    int nneg=params->balancer_num_negative_samples;
    float upper_bounds[NUM_CHANNELS][npos];
    float lower_bounds[NUM_CHANNELS][nneg];
    int channel, rank, cycle;

    for (channel = 0; channel < NUM_CHANNELS; channel++) {
        for (rank = 0; rank < npos; rank++)
            upper_bounds[channel][rank] = -INFINITY;
        for (rank = 0; rank < nneg; rank++)
            lower_bounds[channel][rank] = INFINITY;
    }

    for (cycle = params->balancer_start; cycle < params->balancer_end; cycle++) {
        float signals[NUM_CHANNELS];
        
        decrosstalk_intensity(signals, &intensities[cycle], params->colormatrix);

        printf(" DXTKsignal %d - %.2f %.2f %.2f %.2f\n", cycle - params->balancer_start,
                    signals[0], signals[1], signals[2], signals[3]);

        /* TODO: skip cycles to obey the dark cycles parameters. */
        for (channel = 0; channel < NUM_CHANNELS; channel++) {
            float sig, tmp;
            
            sig = signals[channel];
            for (rank = 0; rank < npos; rank++)
                if (sig > upper_bounds[channel][rank]) {
                    tmp = upper_bounds[channel][rank];
                    upper_bounds[channel][rank] = sig;
                    sig = tmp;
                }

            sig = signals[channel];
            for (rank = 0; rank < nneg; rank++)
                if (sig < lower_bounds[channel][rank]) {
                    tmp = lower_bounds[channel][rank];
                    lower_bounds[channel][rank] = sig;
                    sig = tmp;
                }
        }
    }

    for (channel = 0; channel < NUM_CHANNELS; channel++) {
        float sigsum;

        sigsum = 0.;
        for (rank = 0; rank < npos; rank++)
            sigsum += upper_bounds[channel][rank];
        signal_range_high[channel] = sigsum / npos;

        sigsum = 0.;
        for (rank = 0; rank < nneg; rank++)
            sigsum += lower_bounds[channel][rank];
        signal_range_low[channel] = sigsum / nneg;
    }

    return 0;
}


static int
check_balancer(float *signal_range_high, float *signal_range_low,
               struct IntensitySet *intensities, const char *seq,
               struct PolyARulerParameters *params,
               int *flags)
{
    printf(" BALSEQ = ");
    if (check_balance_minimum(seq, params) < 0) {
        *flags |= PAFLAG_BALANCER_BIASED;
        printf("\n");
        return -1;
    }
    printf("\n");

    if (probe_signal_ranges(signal_range_high, signal_range_low,
                            intensities, params) < 0) {
        *flags |= PAFLAG_BALANCER_SIGNAL_BAD;
        return -1;
    }

    {
        int i;

        for (i = 0; i < 4; i++)
            printf(" BALANCER %d = (%.2f, %.2f)\n", i,
                signal_range_high[i], signal_range_low[i]);
    }

    return 0;
}

static int
process_polya_signal(struct IntensitySet *intensities, const char *seq,
                     int ncycles, const float *signal_range_high,
                     const float *signal_range_low,
                     struct PolyARulerParameters *params)
{
    return 42;
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


int
measure_polya_length(struct CIFData **intensities,
                     const char *sequence_formatted, int ncycles,
                     uint32_t clusterno, int delimiter_end,
                     struct PolyAFinderParameters *finder_params,
                     struct PolyARulerParameters *ruler_params,
                     int *procflags)
{
    float signal_range_high[NUM_CHANNELS], signal_range_low[NUM_CHANNELS];
    int polya_start, polya_end, polya_len;
    uint32_t polya_ret;

    /* Locate the starting position of poly(A) tail if available */
    polya_ret = find_polya(sequence_formatted + delimiter_end,
                           ncycles - delimiter_end, finder_params);
    polya_start = polya_ret >> 16;
    polya_end = polya_ret & 0xffff;
    polya_len = polya_end - polya_start;

    if (polya_start > 0)
        *procflags |= PAFLAG_HAVE_3P_MODIFICATION;

/* XXX set these constants from arguments */
#define POLYA_SIGNAL_PROC_TRIGGER   10
#define READ2_START                 57 /* 0-based, right-excluded */
#define READ2_END                   308 /* 0-based, right-excluded */
#define READ2_LENGTH                (308-57)
/* XXX */
    /* Check balancer region for all spots including non-poly(A)
     * ones. This can be used to suppress the biased filtering of
     * low-quality spots with poly(A)+ tags against poly(A)- tags.
     */
    {
        size_t balancer_length=ruler_params->balancer_end - ruler_params->balancer_start;
        struct IntensitySet spot_intensities[balancer_length];

        fetch_intensity(spot_intensities, intensities, READ2_START,
                        balancer_length, clusterno);

        if (check_balancer(signal_range_high, signal_range_low,
                           spot_intensities, sequence_formatted + READ2_START,
                           ruler_params, procflags) < 0)
            return -1;
    }

    printf(" PAFOUND: %d-%d // %s\n", polya_start, polya_end,
                                      sequence_formatted + delimiter_end);

    /* Process the signals */
    if (polya_len >= POLYA_SIGNAL_PROC_TRIGGER) {
        size_t insert_len=READ2_END - delimiter_end;
        struct IntensitySet spot_intensities[insert_len];
        int polya_len_from_sig;

        fetch_intensity(spot_intensities, intensities, delimiter_end,
                        insert_len, clusterno);

        polya_len_from_sig = process_polya_signal(spot_intensities,
                sequence_formatted + delimiter_end, insert_len,
                signal_range_high, signal_range_low, ruler_params);
        printf(" --> pA (%d)\n", polya_len_from_sig);
    }

    return 0;
}
