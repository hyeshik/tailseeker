/*
 * cdhit-proxy.c
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

#include <stdint.h>
#include <unistd.h>
#include <pthread.h>
#include <assert.h>
#include "tailseq-dedup-approx.h"
#include "cdhit/cdhit-common.h"

#undef MAX_UAA
#define MAX_UAA 4

extern void setaa_to_na();
Options options;

static int
seqdb_add_sequences(SequenceDB &seq_db, struct deduppool *pool)
{
    struct tagcluster *clstr;
    Sequence *seqentry;
    ssize_t clstr_ix;
    int seqno;

    seq_db.Clear();

    seqno = 0;
    for (clstr_ix = pool->tagcluster_left;
            clstr_ix >= 0; clstr_ix = clstr->next) {
        clstr = pool->tagclusters + clstr_ix;
        seqentry = new Sequence();
        if (seqentry == NULL)
            return -1;

        *seqentry = clstr->umi_rep;
	seqentry->dat_length = seqentry->size;
        seqentry->index = seqno++;

        seq_db.sequences.Append(seqentry); 
    }

    return 0;
}

static int
write_sequences(SequenceDB &seq_db, struct deduppool *pool)
{
    struct tagcluster *oldclstr;
    ssize_t *newclstr_ix, oldclstr_ix, next_ix;
    int nclusters=seq_db.rep_seqs.size();
    int seqno, i;

    /* Allocate buffer for indices for new clusters */
    newclstr_ix = (ssize_t *)calloc(nclusters, sizeof(ssize_t));
    if (newclstr_ix == NULL)
        return -1;

    for (i = 0; i < nclusters; i++)
        newclstr_ix[i] = -1;

    /* Reorganize the tagclusters to new clusters. */
    seqno = 0;
    for (oldclstr_ix = pool->tagcluster_left;
            oldclstr_ix >= 0; oldclstr_ix = next_ix) {
        int cdhit_cluster_id;

        oldclstr = pool->tagclusters + oldclstr_ix;
        assert(seqno == seq_db.sequences[seqno]->index);

        next_ix = oldclstr->next;
        cdhit_cluster_id = seq_db.sequences[seqno]->cluster_id;
        if (newclstr_ix[cdhit_cluster_id] < 0) /* The first member */
            newclstr_ix[cdhit_cluster_id] = oldclstr_ix;
        else {
            struct tagcluster *founder;
            founder = pool->tagclusters + newclstr_ix[cdhit_cluster_id];

            /* Replace the representative UMI if new cluster has more clones. */
            if (founder->umi_rep_ndups < oldclstr->umi_rep_ndups) {
                strcpy(founder->umi_rep, oldclstr->umi_rep);
                founder->umi_rep_ndups = oldclstr->umi_rep_ndups;
                memcpy(&founder->trimer_counts.count,
                       &oldclstr->trimer_counts.count, sizeof(oldclstr->trimer_counts));
            }

            /* Concatenate the list of tag alignments to the founder. */
            pool->tagalns[founder->tagaln_tail].next = oldclstr->tagaln_head;
            founder->tagaln_tail = oldclstr->tagaln_tail;
            oldclstr->tagaln_head = oldclstr->tagaln_tail = -1;
            (void)tagcluster_popmiddle(pool, oldclstr);
            tagcluster_free(pool, oldclstr);
        }

        seqno++;
    }

    free((void *)newclstr_ix);
    return 0;
}

int
cdhit_cluster_minitags(struct deduppool *pool)
{
    SequenceDB seq_db;

    seq_db.NAAN = NAAN_array[options.NAA];

    if (seqdb_add_sequences(seq_db, pool) < 0)
        return -1;

    seq_db.SortDivide(options);
    seq_db.DoClustering(options);

    if (write_sequences(seq_db, pool) < 0)
        return -1;

    return 0;
}

int
init_cdhit_internals(double similarity_threshold)
{
    int argc=5;
    const char *argv[argc];

    argv[0] = "tailseq-dedup-approx";
    argv[1] = "-i";
    argv[2] = ".dummy-input";
    argv[3] = "-o";
    argv[4] = ".dummy-output";

    std::cout.setstate(std::ios_base::failbit); // Suppress outputs from cd-hit.

    options.cluster_thd = similarity_threshold;
    options.option_r = 0;
    options.NAA = 10;
    options.NAAN = NAA8;
    options.NAA_top_limit = 12;
    setaa_to_na();
    mat.set_to_na();

    if (options.SetOptions(argc, (char **)argv, false, true) == 0)
        return -1;

    InitNAA(MAX_UAA);

    return 0;
}
