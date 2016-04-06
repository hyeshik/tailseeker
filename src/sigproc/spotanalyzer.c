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


static ssize_t
write_fastq_entry(BGZF *file, const char *entryname, size_t entryname_len,
                  const char *seq, const char *qual, int start, int length)
{
    char writebuf[BUFSIZ];
    char *buf;

    buf = writebuf;

    /* line 1 */
    *buf++ = '@';
    memcpy(buf, entryname, entryname_len);
    buf += entryname_len;
    *buf++ = '\n';

    /* line 2 */
    memcpy(buf, seq + start, length);
    buf += length;
    *buf++ = '\n';

    /* line 3 */
    *buf++ = '+';
    memcpy(buf, entryname, entryname_len);
    buf += entryname_len;
    *buf++ = '\n';

    /* line 4 */
    memcpy(buf, qual + start, length);
    buf += length;
    *buf++ = '\n';

    return bgzf_write(file, writebuf, (size_t)(buf - writebuf));
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
write_taginfo_entry(BGZF *file, struct TailseekerConfig *cfg,
                    struct SampleInfo *sample, uint32_t clusterno,
                    const char *sequence, int delimiter_end,
                    int procflags, int polya_len, int terminal_mods)
{
    char writebuf[BUFSIZ];
    char *buf;
    int i;

    buf = writebuf;

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

    if (BUFSIZ < (int)(buf - writebuf) + sample->umi_total_length) {
        fprintf(stderr, "Tag information buffer is too short.\n");
        return -1;
    }

    for (i = 0; i < sample->umi_ranges_count; i++) {
        struct UMIInterval *umi;

        umi = &sample->umi_ranges[i];
        memcpy(buf, sequence + umi->start, umi->length);
        buf += umi->length;
    }

    *buf++ = '\n';
    *buf = '\0';

    return bgzf_write(file, writebuf, (size_t)(buf - writebuf));
}


static int
write_measurements_to_streams(struct TailseekerConfig *cfg,
                              struct SampleInfo *sample, uint32_t clusterno,
                              const char *sequence_formatted,
                              const char *quality_formatted,
                              int procflags, int delimiter_end, int polya_len,
                              int terminal_mods)
{
#define MAX_ENTRYNAME_LEN   256
    char entryname[MAX_ENTRYNAME_LEN+1];
    size_t entryname_len;

    entryname[MAX_ENTRYNAME_LEN] = '\0';
    snprintf(entryname, MAX_ENTRYNAME_LEN, "%s%04d:%08u:%03d:%03d:",
             cfg->laneid, cfg->tile, (unsigned int)clusterno,
             procflags, polya_len);
    entryname_len = strlen(entryname);

    if (terminal_mods > 0) {
        get_modification_sequence(entryname + entryname_len,
                                  sequence_formatted, delimiter_end,
                                  terminal_mods);
        entryname_len += terminal_mods;
    }

    if (sample->stream_fastq_5 != NULL &&
        write_fastq_entry(sample->stream_fastq_5, entryname,
                          entryname_len,
                          sequence_formatted, quality_formatted,
                          cfg->fivep_start, cfg->fivep_length) < 0)
        return -1;

    if (sample->stream_fastq_3 != NULL) {
        int start, length;

        if (delimiter_end >= 0) {
            start = delimiter_end;
            length = cfg->threep_start + cfg->threep_length - delimiter_end;
        }
        else {
            start = cfg->threep_start;
            length = cfg->threep_length;
        }

        if (write_fastq_entry(sample->stream_fastq_3, entryname,
                              entryname_len,
                              sequence_formatted, quality_formatted,
                              start, length) < 0)
            return -1;
    }

    if (sample->stream_taginfo != NULL &&
        write_taginfo_entry(sample->stream_taginfo, cfg, sample,
                            clusterno, sequence_formatted, delimiter_end,
                            procflags, polya_len, terminal_mods) < 0)
        return -1;

    return 0;
}


int
process_spots(struct TailseekerConfig *cfg, uint32_t firstclusterno,
              struct CIFData **intensities, struct BCLData **basecalls)
{
    uint32_t cycleno, clusterno;
    uint32_t clustersinblock;
    char sequence_formatted[cfg->total_cycles+1], quality_formatted[cfg->total_cycles+1];
    char intensity_formatted[cfg->total_cycles*8+1];
    struct SampleInfo *noncontrol_samples;
    int mismatches;

    int8_t ssw_score_mat[CONTROL_ALIGN_BASE_COUNT * CONTROL_ALIGN_BASE_COUNT];
    int8_t *control_seq;
    ssize_t control_seq_length;
    int32_t min_control_alignment_score, control_alignment_mask_len;


    /* set the starting point of index matching to non-special (other than Unknown and control)
     * samples */
    for (noncontrol_samples = cfg->samples;
         noncontrol_samples != NULL && noncontrol_samples->index[0] != 'X';
         noncontrol_samples = noncontrol_samples->next)
        /* do nothing */;

    /* prepare reference sequence for (PhiX) control */
    control_seq = NULL;
    control_seq_length = -1;
    control_alignment_mask_len = min_control_alignment_score = -1;

    if (cfg->controlinfo.name[0] != '\0') {
        initialize_ssw_score_matrix(ssw_score_mat, CONTROL_ALIGN_MATCH_SCORE,
                                    CONTROL_ALIGN_MISMATCH_SCORE);

        control_seq_length = load_control_sequence(&control_seq);
        if (control_seq_length < 0)
            return -1;

        min_control_alignment_score = cfg->controlinfo.read_length *
                                      CONTROL_ALIGN_MINIMUM_SCORE;
        control_alignment_mask_len = cfg->controlinfo.read_length / 2;
    }

    clustersinblock = intensities[0]->nclusters;
    for (cycleno = 0; cycleno < cfg->total_cycles; cycleno++)
        if (clustersinblock != intensities[cycleno]->nclusters ||
                clustersinblock != basecalls[cycleno]->nclusters) {
            fprintf(stderr, "Inconsistent number of clusters in cycle %d.\n", cycleno);
            return -1;
        }

    mismatches = 0;

    for (clusterno = 0; clusterno < clustersinblock; clusterno++) {
        struct SampleInfo *bc;
        int delimiter_end, procflags=0;
        int polya_len, terminal_mods=-1;

        format_basecalls(sequence_formatted, quality_formatted, basecalls, cfg->total_cycles, clusterno);
        format_intensity(intensity_formatted, intensities, cfg->total_cycles, clusterno, 0);

        bc = assign_barcode(sequence_formatted + cfg->index_start, cfg->index_length,
                            noncontrol_samples, &mismatches);
        if (bc != NULL)
            /* barcode is assigned to a regular sample. do nothing here. */;
        else if (cfg->controlinfo.name[0] == '\0') /* no control sequence is given. treat it Unknown. */
            bc = cfg->samples; /* the first samples in the list is "Unknown". */
        else
            switch (try_alignment_to_control(sequence_formatted, control_seq,
                                             control_seq_length, &cfg->controlinfo,
                                             ssw_score_mat, min_control_alignment_score,
                                             control_alignment_mask_len)) {
                case 0: /* not aligned to control, set as Unknown. */
                    bc = cfg->samples;
                    break;
                case 1: /* aligned. set as control. */
                    bc = cfg->controlinfo.barcode; /* set as control */
                    break;
                case -1: /* error */
                default:
                    fprintf(stderr, "Failed to align read sequence to control.\n");
                    free(control_seq);
                    return -1;
            }

        if (mismatches <= 0) /* no mismatches or falling back to PhiX/Unknown. */
            bc->clusters_mm0++;
        else {
            procflags |= PAFLAG_BARCODE_HAS_MISMATCHES;

            if (mismatches == 1)
                bc->clusters_mm1++;
            else
                bc->clusters_mm2plus++;
        }

        /* Check fingerprint sequences with defined allowed mismatches. */
        if (bc->fingerprint_length > 0) {
            mismatches = count_fingerprint_mismatches(sequence_formatted,
                                                      cfg->threep_start, bc);
            if (mismatches > bc->maximum_fingerprint_mismatches) {
                bc->clusters_fpmismatch++;
                continue;
            }
        }

        if (bc->delimiter_length <= 0)
            polya_len = delimiter_end = -1;
        else {
            delimiter_end = find_delimiter_end_position(sequence_formatted,
                                                        bc, &procflags);
            if (delimiter_end < 0) {
                bc->clusters_nodelim++;
                if (!cfg->keep_no_delimiter)
                    continue;

                polya_len = -1;
            }
            else
                polya_len = measure_polya_length(cfg, intensities,
                                sequence_formatted, clusterno,
                                delimiter_end, &procflags,
                                &terminal_mods);
        }

        if (write_measurements_to_streams(cfg, bc, firstclusterno + clusterno,
                sequence_formatted, quality_formatted,
                procflags, delimiter_end, polya_len, terminal_mods) < 0) {
            perror("process_spots");

            if (control_seq != NULL)
                free(control_seq);
            return -1;
        }
    }

    if (control_seq != NULL)
        free(control_seq);

    return 0;
}

