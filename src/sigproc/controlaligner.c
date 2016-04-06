/*
 * controlaligner.c
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
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <endian.h>
#include <zlib.h>
#include "tailseq-sigproc.h"
#include "contrib/ssw.h"


#define CONTROL_SEQUENCE_SPACING            20 /* space between forward and reverse strands */


static const int8_t DNABASE2NUM[128] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};


static void
initialize_ssw_score_matrix(int8_t *score_mat, int8_t match_score, int8_t mismatch_score)
{
    int l, m, k;

    /* initialize score matrix for Smith-Waterman alignment */
    for (l = k = 0; l < 4; l++) {
        for (m = 0; m < 4; m++)
            score_mat[k++] = (l == m) ? match_score : -mismatch_score;
        score_mat[k++] = 0; /* no penalty for ambiguous base */
    }

    for (m = 0; m < 5; m++)
        score_mat[k++] = 0;
}


static ssize_t
load_control_sequence(int8_t **control_seq)
{
    int8_t *ctlseq, *pctlseq;
    ssize_t len, i;

    len = strlen(phix_control_sequence);
    ctlseq = malloc(len * 2 + CONTROL_SEQUENCE_SPACING);
    if (ctlseq == NULL) {
        perror("load_control_sequence");
        return -1;
    }

    if (len != strlen(phix_control_sequence_rev)) {
        fprintf(stderr, "The lengths of forward and reverse strand for PhiX control sequences "
                        "don't match.\n");
        return -1;
    }

    /* forward strand */
    for (i = 0, pctlseq = ctlseq; i < len; i++)
        *pctlseq++ = DNABASE2NUM[(int)phix_control_sequence[i]];

    for (i = 0; i < CONTROL_SEQUENCE_SPACING; i++)
        *pctlseq++ = 4;

    /* reverse strand */
    for (i = 0; i < len; i++)
        *pctlseq++ = DNABASE2NUM[(int)phix_control_sequence_rev[i]];

    *control_seq = ctlseq;

    return len * 2 + CONTROL_SEQUENCE_SPACING;
}


int
try_alignment_to_control(struct ControlFilterInfo *control_info,
                         const char *sequence_read)
{
    s_profile *alnprof;
    s_align *alnresult;
    int8_t read_seq[control_info->read_length];
    size_t i, j;
    int r;

    /* fast path for perfect matches */
    if (my_strnstr(phix_control_sequence, sequence_read, control_info->read_length) != NULL ||
        my_strnstr(phix_control_sequence_rev, sequence_read, control_info->read_length) != NULL)
        return 1;

    for (i = 0, j = control_info->first_cycle; i < control_info->read_length; i++, j++)
        read_seq[i] = DNABASE2NUM[(int)sequence_read[j]];

    alnprof = ssw_init(read_seq, control_info->read_length, control_info->ssw_score_mat, 5, 0);
    if (alnprof == NULL) {
        perror("try_alignment_to_control");
        return -1;
    }

    alnresult = ssw_align(alnprof, control_info->control_seq,
                          control_info->control_seq_length,
                          CONTROL_ALIGN_GAP_OPEN_SCORE,
                          CONTROL_ALIGN_GAP_EXTENSION_SCORE, 2,
                          control_info->min_control_alignment_score,
                          0, control_info->control_alignment_mask_len);
    r = (alnresult != NULL && alnresult->score1 >= control_info->min_control_alignment_score);

    if (alnresult != NULL)
        align_destroy(alnresult);

    init_destroy(alnprof);

    return r;
}


int
initialize_control_aligner(struct ControlFilterInfo *ctlinfo)
{
    ctlinfo->control_seq = NULL;
    ctlinfo->control_seq_length = -1;
    ctlinfo->control_alignment_mask_len = ctlinfo->min_control_alignment_score = -1;

    if (ctlinfo->name[0] != '\0') {
        initialize_ssw_score_matrix(ctlinfo->ssw_score_mat,
                                    CONTROL_ALIGN_MATCH_SCORE,
                                    CONTROL_ALIGN_MISMATCH_SCORE);

        ctlinfo->control_seq_length = load_control_sequence(&ctlinfo->control_seq);
        if (ctlinfo->control_seq_length < 0)
            return -1;

        ctlinfo->min_control_alignment_score =
                ctlinfo->read_length * CONTROL_ALIGN_MINIMUM_SCORE;
        ctlinfo->control_alignment_mask_len = ctlinfo->read_length / 2;
    }

    return 0;
}

void
free_control_aligner(struct ControlFilterInfo *ctlinfo)
{
    if (ctlinfo->control_seq != NULL) {
        free(ctlinfo->control_seq);
        ctlinfo->control_seq = NULL;
    }
}
