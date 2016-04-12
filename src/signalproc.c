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
#include <ctype.h>
#include "tailseq-sigproc.h"

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
check_balance_minimum(const char *seq, const struct BalancerParameters *params)
{
    short found[256];
    char *base;
    int i;

    memset(found, 0, sizeof(short) * 256);

    for (i = params->start; i < params->end; i++) {
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
                    const struct BalancerParameters *bparams)
{
    int npos=bparams->num_positive_samples;
    int nneg=bparams->num_negative_samples;
    float upper_bounds[NUM_CHANNELS][npos];
    float lower_bounds[NUM_CHANNELS][nneg];
    int channel, rank, cycle;

    for (channel = 0; channel < NUM_CHANNELS; channel++) {
        for (rank = 0; rank < npos; rank++)
            upper_bounds[channel][rank] = -INFINITY;
        for (rank = 0; rank < nneg; rank++)
            lower_bounds[channel][rank] = INFINITY;
    }

    for (cycle = bparams->start; cycle < bparams->end; cycle++) {
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


static int
check_balancer(float *signal_range_low, float *signal_range_bandwidth,
               struct IntensitySet *intensities,
               const float *colormatrix, const char *seq,
               struct BalancerParameters *params, int *flags)
{
#ifdef DEBUG_SIGNAL_PROCESSING
    printf(" BALSEQ = ");
#endif
    if (check_balance_minimum(seq, params) < 0) {
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
                            intensities, colormatrix, params) < 0) {
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


static int
process_polya_signal(struct IntensitySet *intensities, int ncycles,
                     const float *signal_range_low,
                     const float *signal_range_bandwidth,
                     struct PolyARulerParameters *params,
                     int *flags)
{
    int cycle, chan, ndarkcycles;
    int score_cumsum, score_cumsum_max, score_argmax_cycle;
    int downhill_ext_width;
    float downhill_ext_ground;

    ndarkcycles = 0;
    score_cumsum = score_cumsum_max = 0;
    score_argmax_cycle = -1;
    downhill_ext_width = 0;
    downhill_ext_ground = NAN;

    for (cycle = 0; cycle < ncycles; cycle++) {
        float signals[NUM_CHANNELS];
        float normsignals[NUM_CHANNELS];
        float signal_sum;
        float entropy_score, t_intensity_score, composite_score;

        decrosstalk_intensity(signals, &intensities[cycle],
                              params->colormatrix);

        /* Skip assigning scores if spot is dark. */
        signal_sum = 0.f;
        for (chan = 0; chan < NUM_CHANNELS; chan++)
            if (signals[chan] > 0.f)
                signal_sum += signals[chan];

        if (signal_sum < params->dark_cycles_threshold) {
            ndarkcycles++;
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
            continue;
        }

        entropy_score = params->maximum_entropy -
                        shannon_entropy(normsignals);
#if NUM_CHANNELS == 4
#define CHANNEL_T   3
        t_intensity_score = params->t_intensity_score[
                (int)(normsignals[CHANNEL_T] * T_INTENSITY_SCORE_BINS)];
#undef CHANNEL_T
#else
#error Unsupported signal channel count.
#endif

        composite_score = entropy_score * t_intensity_score;
        score_cumsum += -1 + (composite_score >= params->polya_score_threshold) * 2;
        if (score_cumsum > score_cumsum_max) {
            score_cumsum_max = score_cumsum;
            score_argmax_cycle = cycle;

            downhill_ext_ground = entropy_score;
            downhill_ext_width = 0;
        }
        else if (!isnan(downhill_ext_ground)) {
            if (downhill_ext_ground > entropy_score) {
                downhill_ext_width++;
                downhill_ext_ground = entropy_score;
            }
            else
                downhill_ext_ground = NAN;
        }

#ifdef DEBUG_SIGNAL_PROCESSING
        printf("(%03d) %6.2f %6.2f %6.2f %6.2f\t%.2f = %5.2f x %5.2f\t%.2f  %d%s %d\n", cycle,
                signals[0], signals[1], signals[2], signals[3],
                composite_score,
                entropy_score, t_intensity_score, normsignals[3], score_cumsum,
                score_argmax_cycle == cycle ? "*" : "",
                downhill_ext_width);
#endif
    }

#ifdef DEBUG_SIGNAL_PROCESSING
    printf("\n");
#endif

    if (ndarkcycles > 0) {
        *flags |= PAFLAG_DARKCYCLE_EXISTS;
        if (ndarkcycles >= params->max_dark_cycles) {
            *flags |= PAFLAG_DARKCYCLE_OVER_THRESHOLD;
            return -1;
        }
    }

    return score_argmax_cycle + 1 +
           (int)roundf(downhill_ext_width * params->downhill_extension_weight);
}


static int
dump_polya_score(struct IntensitySet *intensities, int ncycles,
                 const float *signal_range_low,
                 const float *signal_range_bandwidth,
                 struct PolyARulerParameters *params,
                 float *scores)
{
    int cycle, chan, ndarkcycles;

    ndarkcycles = 0;

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
            continue;
        }

        entropy_score = params->maximum_entropy -
                        shannon_entropy(normsignals);
#if NUM_CHANNELS == 4
#define CHANNEL_T   3
        t_intensity_score = params->t_intensity_score[
                (int)(normsignals[CHANNEL_T] * T_INTENSITY_SCORE_BINS)];
#undef CHANNEL_T
#else
#error Unsupported signal channel count.
#endif

        scores[cycle] = entropy_score * t_intensity_score;
    }

    if (ndarkcycles >= params->max_dark_cycles)
        return -1;

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
measure_polya_length(struct TailseekerConfig *cfg,
                     struct CIFData **intensities,
                     const char *sequence_formatted, uint32_t clusterno,
                     int delimiter_end, int *procflags,
                     int *terminal_mods, int threep_eff_len)
{
    struct BalancerParameters *balancer_params;
    float signal_range_bandwidth[NUM_CHANNELS];
    float signal_range_low[NUM_CHANNELS];
    int polya_start, polya_end, polya_len;
    uint32_t polya_ret;

    balancer_params = &cfg->balancerparams;

    /* Locate the starting position of poly(A) tail if available */
    polya_ret = find_polya(sequence_formatted + delimiter_end,
                           cfg->threep_start + cfg->threep_length -
                           delimiter_end, &cfg->finderparams);
    polya_start = polya_ret >> 16;
    polya_end = polya_ret & 0xffff;
    polya_len = polya_end - polya_start;
    *terminal_mods = polya_start;

    if (polya_start > 0)
        *procflags |= PAFLAG_HAVE_3P_MODIFICATION;
    if (polya_len == 0)
        *procflags |= PAFLAG_NO_POLYA_DETECTED;

    /* Check balancer region for all spots including non-poly(A)
     * ones. This can be used to suppress the biased filtering of
     * low-quality spots with poly(A)+ tags against poly(A)- tags.
     */
    {
        struct IntensitySet spot_intensities[balancer_params->length];

        fetch_intensity(spot_intensities, intensities, 0,
                        balancer_params->length, clusterno);

        if (check_balancer(signal_range_low, signal_range_bandwidth,
                           spot_intensities, cfg->rulerparams.colormatrix,
                           sequence_formatted + cfg->threep_start,
                           &cfg->balancerparams, procflags) < 0)
            return -1;
    }

#ifdef DEBUG_SIGNAL_PROCESSING
    printf(" PAFOUND: %d-%d // %s\n", polya_start, polya_end,
                                      sequence_formatted + delimiter_end);
#endif

    /* Process the signals */
    if (polya_len >= cfg->finderparams.sigproc_trigger_polya_length) {
        size_t insert_len=cfg->threep_start + threep_eff_len - delimiter_end;
        struct IntensitySet spot_intensities[insert_len];
        int polya_len_from_sig;

        fetch_intensity(spot_intensities, intensities,
                        delimiter_end - cfg->threep_start,
                        insert_len, clusterno);

        polya_len_from_sig = process_polya_signal(spot_intensities,
                insert_len, signal_range_low, signal_range_bandwidth,
                &cfg->rulerparams, procflags);
#ifdef DEBUG_SIGNAL_PROCESSING
        printf(" --> pA (%d)\n", polya_len_from_sig);
#endif

        if (polya_len_from_sig >=
                cfg->finderparams.sigproc_trigger_polya_length) {
            *procflags |= PAFLAG_MEASURED_FROM_FLUORESCENCE;
            return polya_len_from_sig;
        }
        else if (polya_len < cfg->finderparams.naive_ruler_trigger_polya_length)
            return polya_len;
        else {
            /* For ~0.1% long poly(A) tails, fluorescence signal processing
             * gives shorter length than the signal proc trigger. As the
             * dynamic programming-based poly(A) finder gives irrelevantly
             * longer estimation for long poly(A) tails, we use simpler
             * naive ruler.
             */
            const char *ptailend, *ptaildrag, *readend;
            int polya_len_naive, numT;
#define T_COUNTING_WINDOW 2

            ptailend = sequence_formatted + delimiter_end + polya_start;
            readend = sequence_formatted + delimiter_end + insert_len;

            ptailend += cfg->finderparams.naive_ruler_trigger_polya_length;
            polya_len_naive = cfg->finderparams.naive_ruler_trigger_polya_length;
            ptaildrag = ptailend - T_COUNTING_WINDOW;
            numT = T_COUNTING_WINDOW;

            while (ptailend < readend &&
                   (numT + (*ptailend == 'T') >= T_COUNTING_WINDOW)) {
                numT += (*ptailend == 'T') - (*ptaildrag == 'T');
                polya_len_naive += (*ptailend == 'T');

                ptailend++;
                ptaildrag++;
            }

            *procflags |= PAFLAG_MEASURED_USING_NAIVE_RULER;
            return polya_len_naive;
        }
    }

    return polya_len;
}


int
dump_processed_signals(struct TailseekerConfig *cfg, struct SampleInfo *bc,
                       struct CIFData **intensities, const char *sequence_formatted,
                       uint32_t clusterno, int delimiter_end)
/* This function is only used for debugging and optimizing parameters to new platforms,
 * occasionally. Keep the processing algorithm synchronized with `measure_polya_length`.
 */
{
    struct BalancerParameters *balancer_params;
    float signal_range_bandwidth[NUM_CHANNELS];
    float signal_range_low[NUM_CHANNELS];
    int polya_start, polya_end, polya_len, threep_eff_len;
    uint32_t polya_ret;

    balancer_params = &cfg->balancerparams;

    /* Locate the starting position of poly(A) tail if available */
    polya_ret = find_polya(sequence_formatted + delimiter_end,
                           cfg->threep_start + cfg->threep_length -
                           delimiter_end, &cfg->finderparams);
    polya_start = polya_ret >> 16;
    polya_end = polya_ret & 0xffff;
    polya_len = polya_end - polya_start;

    /* Check balancer region for all spots including non-poly(A)
     * ones. This can be used to suppress the biased filtering of
     * low-quality spots with poly(A)+ tags against poly(A)- tags.
     */
    {
        struct IntensitySet spot_intensities[balancer_params->length];
        int dummyflags=0;

        fetch_intensity(spot_intensities, intensities, 0,
                        balancer_params->length, clusterno);

        if (check_balancer(signal_range_low, signal_range_bandwidth,
                           spot_intensities, cfg->rulerparams.colormatrix,
                           sequence_formatted + cfg->threep_start,
                           &cfg->balancerparams, &dummyflags) < 0)
            return 0;
    }

    threep_eff_len = bc->limit_threep_processing;

    /* Process the signals */
    if (polya_len >= cfg->finderparams.sigproc_trigger_polya_length) {
        size_t insert_len=cfg->threep_start + threep_eff_len - delimiter_end;
        struct IntensitySet spot_intensities[insert_len];
        float scores[insert_len];
        ssize_t padding;
        #define ZEROPAD_BLOCK   20
        static const float zeropad[ZEROPAD_BLOCK]={};

        fetch_intensity(spot_intensities, intensities, delimiter_end - cfg->threep_start,
                        insert_len, clusterno);

        if (dump_polya_score(spot_intensities, insert_len, signal_range_low,
                             signal_range_bandwidth, &cfg->rulerparams, scores) < 0)
            return 0;

        pthread_mutex_lock(&bc->statslock);
        if (bgzf_write(bc->stream_signal_dump, scores, sizeof(scores)) < 0) {
            pthread_mutex_unlock(&bc->statslock);
            return -1;
        }

        padding = threep_eff_len - insert_len;
        for (; padding > 0; padding -= ZEROPAD_BLOCK) {
            if (bgzf_write(bc->stream_signal_dump, zeropad,
                    sizeof(float) * (padding < ZEROPAD_BLOCK ? padding : ZEROPAD_BLOCK)) < 0) {
                pthread_mutex_unlock(&bc->statslock);
                return -1;
            }
        }
        pthread_mutex_unlock(&bc->statslock);
    }

    return 0;
}
