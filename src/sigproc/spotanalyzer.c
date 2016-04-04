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


static struct BarcodeInfo *
assign_barcode(const char *indexseq, int barcode_length, struct BarcodeInfo *barcodes,
               int *pmismatches)
{
    struct BarcodeInfo *bestidx, *pidx;
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
find_delimiter_end_position(const char *sequence, struct BarcodeInfo *barcode,
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
count_fingerprint_mismatches(const char *seq, int pos, struct BarcodeInfo *barcodes)
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


int
process_spots(const char *laneid, int tile, int ncycles, uint32_t firstclusterno,
              int threep_start, int threep_length,
              int barcode_start, int barcode_length,
              struct BarcodeInfo *barcodes,
              struct CIFData **intensities, struct BCLData **basecalls,
              struct ControlFilterInfo *control_info,
              struct PolyAFinderParameters *finder_params,
              struct PolyARulerParameters *ruler_params,
              int keep_no_delimiter)
{
    uint32_t cycleno, clusterno;
    uint32_t clustersinblock;
    char sequence_formatted[ncycles+1], quality_formatted[ncycles+1];
    char intensity_formatted[ncycles*8+1];
    struct BarcodeInfo *noncontrol_barcodes;
    int mismatches;

    int8_t ssw_score_mat[CONTROL_ALIGN_BASE_COUNT * CONTROL_ALIGN_BASE_COUNT];
    int8_t *control_seq;
    ssize_t control_seq_length;
    int32_t min_control_alignment_score, control_alignment_mask_len;


    /* set the starting point of index matching to non-special (other than Unknown and control)
     * barcodes */
    for (noncontrol_barcodes = barcodes;
         noncontrol_barcodes != NULL && noncontrol_barcodes->index[0] != 'X';
         noncontrol_barcodes = noncontrol_barcodes->next)
        /* do nothing */;

    /* prepare reference sequence for (PhiX) control */
    control_seq = NULL;
    control_seq_length = -1;
    control_alignment_mask_len = min_control_alignment_score = -1;

    if (control_info != NULL) {
        initialize_ssw_score_matrix(ssw_score_mat, CONTROL_ALIGN_MATCH_SCORE,
                                    CONTROL_ALIGN_MISMATCH_SCORE);

        control_seq_length = load_control_sequence(&control_seq);
        if (control_seq_length < 0)
            return -1;

        min_control_alignment_score = control_info->read_length * CONTROL_ALIGN_MINIMUM_SCORE;
        control_alignment_mask_len = control_info->read_length / 2;
    }

    clustersinblock = intensities[0]->nclusters;
    for (cycleno = 0; cycleno < ncycles; cycleno++)
        if (clustersinblock != intensities[cycleno]->nclusters ||
                clustersinblock != basecalls[cycleno]->nclusters) {
            fprintf(stderr, "Inconsistent number of clusters in cycle %d.\n", cycleno);
            return -1;
        }

    mismatches = 0;

    for (clusterno = 0; clusterno < clustersinblock; clusterno++) {
        struct BarcodeInfo *bc;
        int delimiter_end, procflags=0;
        int polya_len;

        format_basecalls(sequence_formatted, quality_formatted, basecalls, ncycles, clusterno);
        format_intensity(intensity_formatted, intensities, ncycles, clusterno, 0);

        bc = assign_barcode(sequence_formatted + barcode_start, barcode_length,
                            noncontrol_barcodes, &mismatches);
        if (bc != NULL)
            /* barcode is assigned to a regular sample. do nothing here. */;
        else if (control_info == NULL) /* no control sequence is given. treat it Unknown. */
            bc = barcodes; /* the first barcodes in the list is "Unknown". */
        else
            switch (try_alignment_to_control(sequence_formatted, control_seq,
                                             control_seq_length, control_info,
                                             ssw_score_mat, min_control_alignment_score,
                                             control_alignment_mask_len)) {
                case 0: /* not aligned to control, set as Unknown. */
                    bc = barcodes;
                    break;
                case 1: /* aligned. set as control. */
                    bc = control_info->barcode; /* set as control */
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
                                                      threep_start, bc);
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
                if (!keep_no_delimiter)
                    continue;

                polya_len = -1;
            }
            else
                polya_len = measure_polya_length(intensities,
                        sequence_formatted, ncycles, clusterno,
                        threep_start, threep_length,
                        delimiter_end, finder_params, ruler_params,
                        &procflags);
        }

        if (fprintf(bc->stream, "%s%04d\t%d\t%d\t%d\t%d\t",
                    laneid, tile, firstclusterno + clusterno,
                    procflags, delimiter_end, polya_len) < 0) {
            perror("demultiplex_and_write");

            if (control_seq != NULL)
                free(control_seq);
            return -1;
        }
        {
            int i;
            for (i = 0; i < 30; i++)
                fprintf(bc->stream, "%c", sequence_formatted[delimiter_end + i]);
            fprintf(bc->stream, "\n");
        }
/*
        if (fprintf(bc->stream, "%s%04d\t%d\t%d\t%s\t%s\t%s\n", laneid, tile,
                    firstclusterno + clusterno, delimiter_end, sequence_formatted,
                    quality_formatted, intensity_formatted) < 0) {
            perror("demultiplex_and_write");

            if (control_seq != NULL)
                free(control_seq);
            return -1;
        }*/
    }

    if (control_seq != NULL)
        free(control_seq);

    return 0;
}

