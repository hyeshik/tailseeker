/*
 * spotanalyzer.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "tailseq-import.h"

static const char IUPAC_ambiguity_codes[][5] = {
    "A",    /* A */
    "CGT",  /* B */
    "C",    /* C */
    "AGT",  /* D */
    "",     /* E */
    "",     /* F */
    "G",    /* G */
    "ACT",  /* H */
    "",     /* I */
    "",     /* J */
    "GT",   /* K */
    "",     /* L */
    "AC",   /* M */
    "GATC", /* N */
    "",     /* O */
    "",     /* P */
    "",     /* Q */
    "AG",   /* R */
    "CG",   /* S */
    "T",    /* T */
    "T",    /* U */
    "ACG",  /* V */
    "AT",   /* W */
    "GATC", /* X */
    "CT",   /* Y */
    "",     /* Z */
};
#define IUPAC_ambiguity_codes_first     'A'
#define IUPAC_ambiguity_codes_last      'Z'


static struct SampleInfo *
assign_barcode(const char *indexseq, int barcode_length, struct SampleInfo *barcodes,
               int *pmismatches)
{
    struct SampleInfo *bestidx, *pidx;
    int bestmismatches, secondbestfound, i;

    bestidx = NULL;
    bestmismatches = barcode_length + 1;
    secondbestfound = 0;

    /* The first entry in barcodes is "Unknown", thus skip it. */
    for (pidx = barcodes->next; pidx != NULL; pidx = pidx->next) {
        int mismatches=0;

        for (i = 0; i < barcode_length; i++)
            mismatches += (indexseq[i] != pidx->index[i]);

        if (mismatches < bestmismatches) {
            bestidx = pidx;
            bestmismatches = mismatches;
            secondbestfound = 0;
        }
        else if (mismatches == bestmismatches)
            secondbestfound = 1;
    }

    if (bestidx == NULL || secondbestfound ||
            bestmismatches > bestidx->maximum_index_mismatches) {
        *pmismatches = -1;
        return NULL;
    }
    else {
        *pmismatches = bestmismatches;
        return bestidx;
    }
}


static int
find_delimiter_end_position(const char *sequence, struct SampleInfo *barcode,
                            int *flags)
{
    static const int offsets[]={0, -1, 1, 9999};
    const int *poffset;
    const char *delimpos;

    delimpos = sequence + barcode->delimiter_pos;

    for (poffset = &offsets[0]; *poffset < 9999; poffset++) {
        int i, ndiff;
        const char *delimpos_offset;

        delimpos_offset = delimpos + *poffset;

        for (i = ndiff = 0; i < barcode->delimiter_length; i++) {
            const char *matchables;

            if (barcode->delimiter[i] >= IUPAC_ambiguity_codes_first &&
                    barcode->delimiter[i] <= IUPAC_ambiguity_codes_last)
                matchables = IUPAC_ambiguity_codes[
                    barcode->delimiter[i] - IUPAC_ambiguity_codes_first];
            else
                matchables = "";

            for (; *matchables != '\0'; matchables++)
                if (*matchables == delimpos_offset[i])
                    break;

            ndiff += (*matchables == '\0');
        }

        if (ndiff <= barcode->maximum_delimiter_mismatches) {
            if (ndiff > 0)
                *flags |= PAFLAG_DELIMITER_HAS_MISMATCH;
            if (*poffset != 0)
                *flags |= PAFLAG_DELIMITER_IS_SHIFTED;

            return barcode->delimiter_pos + barcode->delimiter_length + *poffset;
        }
    }

    *flags |= PAFLAG_DELIMITER_NOT_FOUND;

    return -1;
}


static int
count_fingerprint_mismatches(const char *seq, int pos, struct SampleInfo *barcodes)
{
    const char *readp, *fpp;
    int mismatches;

    if (barcodes->fingerprint_length <= 0)
        return 0;

    readp = seq + pos;
    fpp = barcodes->fingerprint;
    for (mismatches = 0; *fpp != '\0'; fpp++, readp++)
        mismatches += (*fpp != *readp);

    return mismatches;
}


static int
check_balancer_basecall_quality(struct TailseekerConfig *cfg,
                                struct SampleInfo *sample,
                                const char *phredscore, int *procflags)
{
#define PHRED_BASE  33
    struct BalancerParameters *bparams;
    int qualsum, j;

    bparams = &cfg->balancerparams;

    if (bparams->min_bases_passes > 0) {
        qualsum = 0;

        for (j = bparams->start; j < bparams->end; j++)
            qualsum += ((phredscore[j] - PHRED_BASE) >= bparams->min_quality);

        if (qualsum < bparams->min_bases_passes) {
            *procflags |= PAFLAG_BALANCER_CALL_QUALITY_BAD;

            pthread_mutex_lock(&sample->statslock);
            sample->clusters_qcfailed++;
            pthread_mutex_unlock(&sample->statslock);

            if (!cfg->keep_low_quality_balancer)
                return -1;
        }
    }

    return 0;
}


static ssize_t
write_seqqual_entry(char **pbuffer, uint32_t clusterno,
                    const char *seq, const char *qual,
                    int start_5p, int length_5p, int start_3p, int length_3p)
{
    ssize_t written;
    char *buf;

    buf = *pbuffer;

    /* clusterno */
    written = sprintf(buf, "%u\t", (unsigned int)clusterno);
    buf += written;

    /* 5'-side read sequence */
    memcpy(buf, seq + start_5p, length_5p);
    buf += length_5p;
    *buf++ = '\t';

    /* 5'-side read quality */
    memcpy(buf, qual + start_5p, length_5p);
    buf += length_5p;
    *buf++ = '\t';

    /* 3'-side read sequence */
    memcpy(buf, seq + start_3p, length_3p);
    buf += length_3p;
    *buf++ = '\t';

    /* 3'-side read quality */
    memcpy(buf, qual + start_3p, length_3p);
    buf += length_3p;
    *buf++ = '\n';

    written = buf - *pbuffer;
    *pbuffer = buf;

    return written;
}


static void
get_modification_sequence(char *dest, const char *seq, int delimiter_end,
                          int terminal_mods)
{
    if (terminal_mods >= 1) {
        const char *seqptr;
        int i;

        seqptr = seq + delimiter_end + terminal_mods - 1;
        for (i = 0; i < terminal_mods; i++, seqptr--)
            switch (*seqptr) {
            case 'A': *dest++ = 'T'; break;
            case 'C': *dest++ = 'G'; break;
            case 'G': *dest++ = 'C'; break;
            case 'T': *dest++ = 'A'; break;
            default:  *dest++ = 'N'; break;
            }
    }

    *dest = '\0';
}


static int
write_taginfo_entry(char **pbuffer, struct TailseekerConfig *cfg,
                    struct SampleInfo *sample, uint32_t clusterno,
                    const char *sequence, int delimiter_end,
                    int procflags, int polya_len, int terminal_mods)
{
    ssize_t written;
    char *buf;
    int i;

    buf = *pbuffer;

    sprintf(buf, "%u\t%d\t%d\t", (unsigned int)clusterno, procflags,
            polya_len);
    buf += strlen(buf);

    if (terminal_mods > 0) {
        get_modification_sequence(buf, sequence, delimiter_end,
                                  terminal_mods);
        buf += terminal_mods;
    }
    *buf++ = '\t';

    for (i = 0; i < sample->umi_ranges_count; i++) {
        struct UMIInterval *umi;

        umi = &sample->umi_ranges[i];
        memcpy(buf, sequence + umi->start, umi->length);
        buf += umi->length;
    }

    *buf++ = '\n';

    written = buf - *pbuffer;
    *pbuffer = buf;

    return written;
}


static int
write_measurements_to_buffers(struct TailseekerConfig *cfg,
                              struct WriteBuffer *wbuf,
                              struct SampleInfo *sample, uint32_t clusterno,
                              const char *sequence_formatted,
                              const char *quality_formatted,
                              int procflags, int delimiter_end, int polya_len,
                              int terminal_mods)
{
    struct WriteBuffer *wb;

    wb = &wbuf[sample->numindex];

    if (sample->stream_seqqual != NULL) {
        int start_3p, length_3p;

        if (delimiter_end >= 0) {
            start_3p = delimiter_end;
            length_3p = cfg->threep_start + cfg->threep_length - delimiter_end;
        }
        else {
            start_3p = cfg->threep_start;
            length_3p = cfg->threep_length;
        }

        if (length_3p > cfg->threep_seqqual_output_length)
            length_3p = cfg->threep_seqqual_output_length;

        if (write_seqqual_entry(&wb->buf_seqqual, clusterno,
                                sequence_formatted, quality_formatted,
                                cfg->fivep_start, cfg->fivep_length,
                                start_3p, length_3p) < 0)
            return -1;
    }

    if (sample->stream_taginfo != NULL &&
        write_taginfo_entry(&wb->buf_taginfo, cfg, sample,
                            clusterno, sequence_formatted, delimiter_end,
                            procflags, polya_len, terminal_mods) < 0)
        return -1;

    return 0;
}


static ssize_t
sync_write_out_buffer(BGZF *stream, const char *content, size_t size,
                      struct WriteHandleSync *sync, int jobid)
{
    ssize_t written;

    while (1) {
        pthread_mutex_lock(&sync->lock);
        if (sync->jobs_written >= jobid) {
            pthread_mutex_unlock(&sync->lock);
            break;
        }

        pthread_cond_wait(&sync->wakeup, &sync->lock);
        pthread_mutex_unlock(&sync->lock);
    }

    if (size > 0) {
        written = bgzf_write(stream, content, size);
        if (written < 0)
            return -1;
    }
    else
        written = 0;

    pthread_mutex_lock(&sync->lock);
    sync->jobs_written++;
    pthread_cond_broadcast(&sync->wakeup);
    pthread_mutex_unlock(&sync->lock);

    return written;
}


static int
write_polya_score(struct SampleInfo *sample, const float *score, int length,
                  const char *downhill, uint32_t clusterno, int first_cycle)
{
    static const signal_packet_t zeropad[256] = {{0, 0}};
    struct SignalRecordHeader header;
    signal_packet_t sigscores[length];
    unsigned int s;
    int i;

    if (length > sample->signal_dump_length)
        length = sample->signal_dump_length;

    header.clusterno = clusterno;
    header.first_cycle = first_cycle;
    header.valid_cycle_count = length;

    memset(sigscores, 0, sizeof(signal_packet_t) * length);
    for (i = 0; i < length; i++)
        if (!isnan(score[i])) {
            assert(score[i] >= 0.f && score[i] <= 1.f);
            s = 1 + (int)(score[i] * (float)(SIGNALPACKET_SCORE_MAX - 1));
            if (s >= SIGNALPACKET_SCORE_MAX)
                s = SIGNALPACKET_SCORE_MAX;
            sigscores[i].score = s;
            sigscores[i].downhill = downhill[i];
        }

    pthread_mutex_lock(&sample->signal_writer_lock);

    if (bgzf_write(sample->stream_signal, (void *)&header,
                sizeof(header)) < 0 ||
            bgzf_write(sample->stream_signal, (void *)sigscores,
                    sizeof(signal_packet_t) * length) < 0) {
        perror("write_polya_score");
        return -1;
    }

    if (length < sample->signal_dump_length &&
            bgzf_write(sample->stream_signal, (void *)zeropad,
                       sizeof(signal_packet_t) *
                       (sample->signal_dump_length - length)) < 0) {
        perror("write_polya_score");
        return -1;
    }

    pthread_mutex_unlock(&sample->signal_writer_lock);

    return 0;
}


static void
add_polya_score_sample(cluster_count_t *samplecounts, float *scores,
                       int length, int firstcycle, int sampling_bins)
{
    cluster_count_t *sptr;
    int i, binno;

    sptr = samplecounts + (firstcycle * sampling_bins);
    for (i = 0; i < length; i++, scores++, sptr += sampling_bins)
        if (isfinite(*scores)) {
            binno = *scores * sampling_bins;
            if (binno >= sampling_bins)
                binno = sampling_bins - 1;

            sptr[binno]++;
        }
}


static int
check_sequence_for_fair_sampling(struct FairSamplingCount *fair_sampling,
                                 struct PolyASeederParameters *params,
                                 const char *sequence)
{
    const char *end, *hashptr;
    uint64_t v;
    int r;

    end = sequence + params->fair_sampling_fingerprint_length;

    v = 0;
    for (hashptr = sequence; hashptr < end; hashptr++)
        v = ((v << 3) | ((*hashptr & 7) ^ (v >> 60))) & ((1UL << 63) - 1UL);
    v %= params->fair_sampling_hash_space_size;

    pthread_mutex_lock(&fair_sampling->lock);

    if (params->fair_sampling_max_count == 0 ||
            fair_sampling->count[v] < params->fair_sampling_max_count) {
        fair_sampling->count[v]++;
        r = 0;
    }
    else
        r = -1;

    pthread_mutex_unlock(&fair_sampling->lock);

    return r;
}


static int
process_polya_signal(struct TailseekerConfig *cfg, uint32_t clusterno,
                     uint32_t global_clusterno, struct SampleInfo *sample,
                     const char *sequence,
                     struct CIFData **intensities, int delimiter_end,
                     int *terminal_mods,
                     cluster_count_t *pos_score_counts,
                     cluster_count_t *neg_score_counts,
                     struct FairSamplingCount *fair_sampling,
                     int *procflags)
{
    int polya_start, polya_end, balancer_len, insert_len, polya_len;
    uint32_t polya_ret;
    float signal_range_bandwidth[NUM_CHANNELS];
    float signal_range_low[NUM_CHANNELS];
    struct PolyASeederParameters params;

    memcpy(&params, &cfg->seederparams, sizeof(params));

    /* Locate the starting position of poly(A) tail if available */
    polya_ret = find_polya(sequence + delimiter_end,
                           cfg->threep_start + cfg->threep_length - delimiter_end,
                           &cfg->finderparams);
    polya_start = polya_ret >> 16;
    polya_end = polya_ret & 0xffff;
    polya_len = polya_end - polya_start;
    *terminal_mods = polya_start;

    balancer_len = delimiter_end - cfg->threep_start;
    if (balancer_len > cfg->balancerparams.length)
        balancer_len = cfg->balancerparams.length;

    if (polya_start > 0)
        *procflags |= PAFLAG_HAVE_3P_MODIFICATION;
    if (polya_len > 0)
        *procflags |= PAFLAG_POLYA_DETECTED;

    /* Check balancer region for all spots including non-poly(A)
     * ones. This can be used to suppress the biased filtering of
     * low-quality spots with poly(A)+ tags against poly(A)- tags.
     */
    {
        struct IntensitySet spot_intensities[balancer_len];

        fetch_intensity(spot_intensities, intensities, 0, balancer_len, clusterno);

        if (check_balancer(signal_range_low, signal_range_bandwidth,
                           spot_intensities, cfg->rulerparams.colormatrix,
                           sequence + cfg->threep_start,
                           &cfg->balancerparams, balancer_len, procflags) < 0)
            return -1;
    }

    insert_len = cfg->threep_start + cfg->threep_length - delimiter_end;
    
    {
        struct IntensitySet spot_intensities[insert_len];
        float scores[insert_len], contrast_score;
        char downhill[insert_len];
        int max_contrast_pos, scan_len;

        /* Fetch signal intensities and calculate poly(A) scores. */
        fetch_intensity(spot_intensities, intensities,
                        delimiter_end - cfg->threep_start,
                        insert_len, clusterno);

        if (compute_polya_score(spot_intensities, insert_len,
                                signal_range_low, signal_range_bandwidth,
                                &cfg->rulerparams, scores, downhill,
                                procflags) < 0)
            return -1;

        scan_len = insert_len - polya_start;
        if (scan_len <= 0)
            return polya_len;

        /* Evaluate the signal trends if it stems from a poly(A) tail.
         * A poly(A) tail has a great contrast among score values
         * divided by a certain point. */
        if (polya_len >= params.seed_trigger_polya_length) {
            max_contrast_pos = find_max_cumulative_contrast(scores + polya_start,
                                    scan_len, params.max_cctr_scan_left_space,
                                    params.max_cctr_scan_right_space,
                                    &contrast_score);

            if (max_contrast_pos >= params.polya_boundary_pos &&
                    contrast_score >= params.required_cdf_contrast) {
                int sampling_len = -polya_start;
                sampling_len += min_int(scan_len,
                        params.polya_boundary_pos - params.polya_sampling_gap);
                if (sampling_len > 0)
                    add_polya_score_sample(pos_score_counts,
                                           scores + polya_start,
                                           sampling_len,
                                           delimiter_end + polya_start,
                                           params.dist_sampling_bins);
            }
        }
        else if (polya_len <= params.negative_sample_polya_length &&
                    scan_len >= params.fair_sampling_fingerprint_length &&
                    check_sequence_for_fair_sampling(fair_sampling,
                        &params, sequence + delimiter_end + polya_start) == 0)
            add_polya_score_sample(neg_score_counts, scores + polya_start,
                                   scan_len, delimiter_end + polya_start,
                                   params.dist_sampling_bins);

        /* Write computed poly(A) scores of long poly(A) candidates
         * for later evaluation. */
        if (polya_len >= cfg->finderparams.sigproc_trigger_polya_length) {
            if (write_polya_score(sample, scores + polya_start, polya_len,
                                  downhill, global_clusterno,
                                  delimiter_end + polya_start) < 0) {
                fprintf(stderr, "Failed to write a poly(A) score.\n");
                return polya_len; //Modified by hanju: apply polya_len instead of scan_len
            }
        }
        else if (polya_len == 0) {
            if (write_polya_score(sample, scores + polya_start, scan_len,
                                  downhill, global_clusterno,
                                  delimiter_end + polya_start) < 0) {
                fprintf(stderr, "Failed to write a poly(A) score.\n");
                return polya_len; //Modified by hanju: write polya_len even undetectable
            }
        }
    }

    return polya_len;
}


int
process_spots(struct TailseekerConfig *cfg, uint32_t firstclusterno,
              struct CIFData **intensities, struct BCLData **basecalls,
              struct WriteBuffer *wbuf0, struct WriteBuffer *wbuf,
              cluster_count_t *pos_score_counts,
              cluster_count_t *neg_score_counts,
              struct FairSamplingCount *fair_sampling,
              int jobid, uint32_t cln_start, uint32_t cln_end)
{
    uint32_t clusterno;
    char sequence_formatted[cfg->total_cycles+1], quality_formatted[cfg->total_cycles+1];
    struct SampleInfo *noncontrol_samples;
    int mismatches;

    /* set the starting point of index matching to non-special (other than Unknown and control)
     * samples */
    for (noncontrol_samples = cfg->samples;
         noncontrol_samples != NULL && noncontrol_samples->index[0] != 'X';
         noncontrol_samples = noncontrol_samples->next)
        /* do nothing */;

    mismatches = 0;

    for (clusterno = cln_start; clusterno < cln_end; clusterno++) {
        struct SampleInfo *sample;
        int delimiter_end, procflags=0;
        int polya_status, terminal_mods=-1;

        format_basecalls(sequence_formatted, quality_formatted, basecalls,
                         cfg->total_cycles, clusterno);

        sample = assign_barcode(sequence_formatted + cfg->index_start, cfg->index_length,
                                noncontrol_samples, &mismatches);
        if (sample != NULL)
            /* barcode is assigned to a regular sample. do nothing here. */;
        else if (cfg->controlinfo.name[0] == '\0') /* no control sequence is given. treat it Unknown. */
            sample = cfg->samples; /* the first samples in the list is "Unknown". */
        else
            switch (try_alignment_to_control(&cfg->controlinfo, sequence_formatted)) {
                case 0: /* not aligned to control, set as Unknown. */
                    sample = cfg->samples;
                    break;
                case 1: /* aligned. set as control. */
                    sample = cfg->controlinfo.barcode; /* set as control */
                    break;
                case -1: /* error */
                default:
                    fprintf(stderr, "Failed to align read sequence to control.\n");
                    return -1;
            }

        pthread_mutex_lock(&sample->statslock);
        if (mismatches <= 0) /* no mismatches or falling back to PhiX/Unknown. */
            sample->clusters_mm0++;
        else {
            procflags |= PAFLAG_BARCODE_HAS_MISMATCHES;

            if (mismatches == 1)
                sample->clusters_mm1++;
            else
                sample->clusters_mm2plus++;
        }
        pthread_mutex_unlock(&sample->statslock);

        /* Check fingerprint sequences with defined allowed mismatches. */
        if (sample->fingerprint_length > 0) {
            mismatches = count_fingerprint_mismatches(sequence_formatted,
                                                      sample->fingerprint_pos, sample);
            if (mismatches > sample->maximum_fingerprint_mismatches) {
                pthread_mutex_lock(&sample->statslock);
                sample->clusters_fpmismatch++;
                pthread_mutex_unlock(&sample->statslock);
                continue;
            }
        }

        /* Check basecalling quality scores in the balancer in 3'-side read.
         * This will represent how good the signal quality is. Using any
         * among other regions leads to a biased sampling against long poly(A)
         * tails. */
        if (sample->umi_ranges_count > 0 &&
                check_balancer_basecall_quality(cfg, sample, quality_formatted,
                                                &procflags) < 0)
            continue;

        if (sample->delimiter_length <= 0)
            polya_status = delimiter_end = -1;
        else {
            delimiter_end = find_delimiter_end_position(sequence_formatted,
                                                        sample, &procflags);
            if (delimiter_end < 0) {
                pthread_mutex_lock(&sample->statslock);
                sample->clusters_nodelim++;
                pthread_mutex_unlock(&sample->statslock);
                if (!cfg->keep_no_delimiter)
                    continue;

                polya_status = -1;
            }
            else
                polya_status = process_polya_signal(cfg, clusterno,
                      firstclusterno + clusterno, sample,
                      sequence_formatted,
                      intensities, delimiter_end, &terminal_mods,
                      pos_score_counts, neg_score_counts,
                      fair_sampling, &procflags);
        }

        if (write_measurements_to_buffers(cfg, wbuf, sample,
                firstclusterno + clusterno, sequence_formatted, quality_formatted,
                procflags, delimiter_end, polya_status, terminal_mods) < 0) {
            perror("process_spots");

            return -1;
        }
    }

    {
        struct SampleInfo *sample;

        for (sample = cfg->samples; sample != NULL; sample = sample->next) {
            if (sync_write_out_buffer(sample->stream_seqqual,
                                      wbuf0[sample->numindex].buf_seqqual,
                                      (size_t)(wbuf[sample->numindex].buf_seqqual -
                                               wbuf0[sample->numindex].buf_seqqual),
                                      &sample->wsync_seqqual, jobid) < 0)
                return -1;

            if (sync_write_out_buffer(sample->stream_taginfo,
                                      wbuf0[sample->numindex].buf_taginfo,
                                      (size_t)(wbuf[sample->numindex].buf_taginfo -
                                               wbuf0[sample->numindex].buf_taginfo),
                                      &sample->wsync_taginfo, jobid) < 0)
                return -1;
        }
    }

    return 0;
}
