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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include "tailseq-import.h"

//#define DEBUG_SIGNAL_PROCESSING

static inline float
shannon_entropy(const float *probs)
{
    float res;
    int i;

    res = 0.f;
    for (i = 0; i < NUM_CHANNELS; i++)
        if (probs[i] > 0.f)
            res -= probs[i] * logf(probs[i]);

    return res;
}


static float
calculate_maximum_entropy(void)
{
    float prob[NUM_CHANNELS];
    int i;

    for (i = 0; i < NUM_CHANNELS; i++)
        prob[i] = 1. / NUM_CHANNELS;

    return shannon_entropy(prob);
}


void
precalc_score_tables(struct PolyARulerParameters *params, float k, float center)
{
    float vpoint;
    int i;

    params->maximum_entropy = calculate_maximum_entropy();

    for (i = 0; i < T_INTENSITY_SCORE_BINS; i++) {
        vpoint = (1. * (.5 + i)) / T_INTENSITY_SCORE_BINS;
        params->t_intensity_score[i] = 1. /
                    (1. + expf(-k * (vpoint - center)));
    }

    params->t_intensity_score[T_INTENSITY_SCORE_BINS] = 1.f; /* for t=1.0 */
}


static inline void
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
check_balance_minimum(const char *seq, const struct BalancerParameters *params,
                      int balancer_len)
{
    short found[256];
    char *base;
    int i, end;

    memset(found, 0, sizeof(short) * 256);
    end = params->start + balancer_len;

    for (i = params->start; i < end; i++) {
#ifdef DEBUG_SIGNAL_PROCESSING
        printf("%c", seq[i]);
#endif
        found[(int)seq[i]]++;
    }

    for (base = "ACGT"; *base != '\0'; base++)
        if (found[(int)*base] < params->minimum_occurrence)
            return -1;

    return 0;
}


static int
probe_signal_ranges(float *signal_range_low, float *signal_range_bandwidth,
                    const struct IntensitySet *intensities,
                    const float *colormatrix,
                    const struct BalancerParameters *bparams,
                    int balancer_len)
{
    int npos=bparams->num_positive_samples;
    int nneg=bparams->num_negative_samples;
    float upper_bounds[NUM_CHANNELS][npos];
    float lower_bounds[NUM_CHANNELS][nneg];
    int channel, rank, cycle, end;

    for (channel = 0; channel < NUM_CHANNELS; channel++) {
        for (rank = 0; rank < npos; rank++)
            upper_bounds[channel][rank] = -INFINITY;
        for (rank = 0; rank < nneg; rank++)
            lower_bounds[channel][rank] = INFINITY;
    }

    end = bparams->start + balancer_len;
    for (cycle = bparams->start; cycle < end; cycle++) {
        float signals[NUM_CHANNELS];

        decrosstalk_intensity(signals, &intensities[cycle], colormatrix);

#ifdef DEBUG_SIGNAL_PROCESSING
        printf(" DXTKsignal %d - %.2f %.2f %.2f %.2f\n", cycle - bparams->start,
                    signals[0], signals[1], signals[2], signals[3]);
#endif

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

        sigsum = 0.f;
        for (rank = 0; rank < nneg; rank++)
            sigsum += lower_bounds[channel][rank];
        signal_range_low[channel] = sigsum / nneg;

        sigsum = 0.f;
        for (rank = 0; rank < npos; rank++)
            sigsum += upper_bounds[channel][rank];
        signal_range_bandwidth[channel] = sigsum / npos
                                    - signal_range_low[channel];
    }

    return 0;
}


int
check_balancer(float *signal_range_low, float *signal_range_bandwidth,
               struct IntensitySet *intensities,
               const float *colormatrix, const char *seq,
               struct BalancerParameters *params, int balancer_len,
               int *flags)
{
#ifdef DEBUG_SIGNAL_PROCESSING
    printf(" BALSEQ = ");
#endif
    if (check_balance_minimum(seq, params, balancer_len) < 0) {
        *flags |= PAFLAG_BALANCER_BIASED;
#ifdef DEBUG_SIGNAL_PROCESSING
        printf("\n");
#endif
        return -1;
    }
#ifdef DEBUG_SIGNAL_PROCESSING
    printf("\n");
#endif

    if (probe_signal_ranges(signal_range_low, signal_range_bandwidth,
                            intensities, colormatrix, params,
                            balancer_len) < 0) {
        *flags |= PAFLAG_BALANCER_SIGNAL_BAD;
        return -1;
    }

#ifdef DEBUG_SIGNAL_PROCESSING
    {
        int i;

        for (i = 0; i < 4; i++)
            printf(" BALANCER %d = (%.2f, %.2f)\n", i,
                signal_range_low[i], signal_range_bandwidth[i]);
    }
#endif

    return 0;
}


static inline void
normalize_signals(float *normalized, const float *signals,
                  const float *signal_range_low,
                  const float *signal_range_bandwidth)
{
    float signal_sum;
    int chan;

    signal_sum = 0.f;
    for (chan = 0; chan < NUM_CHANNELS; chan++) {
        normalized[chan] = (signals[chan] - signal_range_low[chan])
                            / signal_range_bandwidth[chan];
        if (normalized[chan] < 0.f)
            normalized[chan] = 0.f;

        signal_sum += normalized[chan];
    }

    for (chan = 0; chan < NUM_CHANNELS; chan++)
        normalized[chan] /= signal_sum;
}


int
compute_polya_score(struct IntensitySet *intensities, int ncycles,
                    const float *signal_range_low,
                    const float *signal_range_bandwidth,
                    struct PolyARulerParameters *params,
                    float *scores, char *downhill, int *procflags)
{
    int cycle, chan, ndarkcycles;
    float entropy_prev;

    ndarkcycles = 0;
    entropy_prev = NAN;

    for (cycle = 0; cycle < ncycles; cycle++) {
        float signals[NUM_CHANNELS];
        float normsignals[NUM_CHANNELS];
        float signal_sum;
        float entropy_score, t_intensity_score;

        decrosstalk_intensity(signals, &intensities[cycle],
                              params->colormatrix);

        /* Skip assigning scores if spot is dark. */
        signal_sum = 0.f;
        for (chan = 0; chan < NUM_CHANNELS; chan++)
            if (signals[chan] > 0.f)
                signal_sum += signals[chan];

        if (signal_sum < params->dark_cycles_threshold) {
            ndarkcycles++;
            *scores++ = entropy_prev = NAN;
            *downhill++ = 0;
            continue;
        }

        /* Adjust signals to fit in the spot dynamic range */
        normalize_signals(normsignals, signals, signal_range_low,
                          signal_range_bandwidth);
        if (isnan(normsignals[0])) {
            /* Normalized signals can be NaNs altogether if intensities of
             * all channels are zero or negative after normalization.
             * It is treated as a dark cycle in this case. */
            ndarkcycles++;
            *scores++ = entropy_prev = NAN;
            *downhill++ = 0;
            continue;
        }

        entropy_score = 1.f - shannon_entropy(normsignals) /
                        params->maximum_entropy;
#if NUM_CHANNELS == 4
#define CHANNEL_T   3
        t_intensity_score = params->t_intensity_score[
                (int)(normsignals[CHANNEL_T] * T_INTENSITY_SCORE_BINS)];
#undef CHANNEL_T
#else
#error Unsupported signal channel count.
#endif

        *scores++ = entropy_score * t_intensity_score;
        *downhill++ = (char)(entropy_score < entropy_prev);
        entropy_prev = entropy_score;
    }

    if (ndarkcycles > 0) {
        *procflags |= PAFLAG_DARKCYCLE_EXISTS;
        if (ndarkcycles >= params->max_dark_cycles) {
            *procflags |= PAFLAG_DARKCYCLE_OVER_THRESHOLD;
            return -1;
        }
    }

    return 0;
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
find_max_cumulative_contrast(const float *scores, int length, int leftspace,
                             int rightspace, float *pmax_score)
{
    double score_cumsum[length], score_total;
    double max_score, score;
    int i, max_contrast_boundary, right;

    score_total = 0.;
    for (i = 0; i < length; i++) {
        score_total += scores[i];
        score_cumsum[i] = score_total;
    }

    assert(leftspace >= 1 && rightspace >= 1);
    if (leftspace + rightspace >= length)
        return -1;

    max_contrast_boundary = leftspace - 1;
    max_score = score_cumsum[max_contrast_boundary] / max_contrast_boundary -
                (score_total - score_cumsum[max_contrast_boundary]) /
                (length - max_contrast_boundary);

    right = length - rightspace;
    for (i = leftspace; i <= right; i++) {
        score = score_cumsum[i] / i - (score_total - score_cumsum[i]) / (length - i);
        if (score > max_score) {
            max_contrast_boundary = i;
            max_score = score;
        }
    }

    *pmax_score = max_score;
    return max_contrast_boundary + 1;
}
