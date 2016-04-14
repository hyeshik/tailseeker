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

#define _BSD_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "tailseq-sigproc.h"


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

        for (i = ndiff = 0; i < barcode->delimiter_length; i++)
            ndiff += (barcode->delimiter[i] != delimpos_offset[i]);

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
write_seqqual_entry(char **pbuffer, struct TailseekerConfig *cfg, uint32_t clusterno,
                    const char *seq, const char *qual,
                    int start_5p, int length_5p, int start_3p, int length_3p)
{
    ssize_t written;
    char *buf;

    buf = *pbuffer;

    /* lane id and clusterno */
    written = sprintf(buf, "%s%04d\t%u\t", cfg->laneid, cfg->tile, (unsigned int)clusterno);
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

    sprintf(buf, "%s%04d\t%u\t%d\t%d\t",
             cfg->laneid, cfg->tile, (unsigned int)clusterno,
             procflags, polya_len);
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

        if (length_3p > cfg->threep_output_length)
            length_3p = cfg->threep_output_length;

        if (write_seqqual_entry(&wb->buf_seqqual, cfg, clusterno,
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


int
process_spots(struct TailseekerConfig *cfg, uint32_t firstclusterno,
              struct CIFData **intensities, struct BCLData **basecalls,
              struct WriteBuffer *wbuf0, struct WriteBuffer *wbuf,
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
        struct SampleInfo *bc;
        int delimiter_end, procflags=0;
        int polya_len, terminal_mods=-1;

        format_basecalls(sequence_formatted, quality_formatted, basecalls,
                         cfg->total_cycles, clusterno);

        bc = assign_barcode(sequence_formatted + cfg->index_start, cfg->index_length,
                            noncontrol_samples, &mismatches);
        if (bc != NULL)
            /* barcode is assigned to a regular sample. do nothing here. */;
        else if (cfg->controlinfo.name[0] == '\0') /* no control sequence is given. treat it Unknown. */
            bc = cfg->samples; /* the first samples in the list is "Unknown". */
        else
            switch (try_alignment_to_control(&cfg->controlinfo, sequence_formatted)) {
                case 0: /* not aligned to control, set as Unknown. */
                    bc = cfg->samples;
                    break;
                case 1: /* aligned. set as control. */
                    bc = cfg->controlinfo.barcode; /* set as control */
                    break;
                case -1: /* error */
                default:
                    fprintf(stderr, "Failed to align read sequence to control.\n");
                    return -1;
            }

        pthread_mutex_lock(&bc->statslock);
        if (mismatches <= 0) /* no mismatches or falling back to PhiX/Unknown. */
            bc->clusters_mm0++;
        else {
            procflags |= PAFLAG_BARCODE_HAS_MISMATCHES;

            if (mismatches == 1)
                bc->clusters_mm1++;
            else
                bc->clusters_mm2plus++;
        }
        pthread_mutex_unlock(&bc->statslock);

        /* Check fingerprint sequences with defined allowed mismatches. */
        if (bc->fingerprint_length > 0) {
            mismatches = count_fingerprint_mismatches(sequence_formatted,
                                                      bc->fingerprint_pos, bc);
            if (mismatches > bc->maximum_fingerprint_mismatches) {
                pthread_mutex_lock(&bc->statslock);
                bc->clusters_fpmismatch++;
                pthread_mutex_unlock(&bc->statslock);
                continue;
            }
        }

        /* Check basecalling quality scores in the balancer in 3'-side read.
         * This will represent how good the signal quality is. Using any
         * among other regions leads to a biased sampling against long poly(A)
         * tails. */
        if (bc->umi_ranges_count > 0 &&
                check_balancer_basecall_quality(cfg, bc, quality_formatted,
                                                &procflags) < 0)
            continue;

        if (bc->delimiter_length <= 0)
            polya_len = delimiter_end = -1;
        else {
            delimiter_end = find_delimiter_end_position(sequence_formatted,
                                                        bc, &procflags);
            if (delimiter_end < 0) {
                pthread_mutex_lock(&bc->statslock);
                bc->clusters_nodelim++;
                pthread_mutex_unlock(&bc->statslock);
                if (!cfg->keep_no_delimiter)
                    continue;

                polya_len = -1;
            }
            else
                polya_len = measure_polya_length(cfg, intensities,
                                sequence_formatted, clusterno,
                                delimiter_end, &procflags,
                                &terminal_mods,
                                bc->limit_threep_processing);
        }

        if (write_measurements_to_buffers(cfg, wbuf, bc,
                firstclusterno + clusterno, sequence_formatted, quality_formatted,
                procflags, delimiter_end, polya_len, terminal_mods) < 0) {
            perror("process_spots");

            return -1;
        }

        if (bc->stream_signal_dump != NULL &&
                dump_processed_signals(cfg, bc, intensities, sequence_formatted,
                                       clusterno, delimiter_end) < 0) {
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
