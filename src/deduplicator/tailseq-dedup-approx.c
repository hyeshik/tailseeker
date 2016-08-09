/*
 * tailseq-dedup-approx.c
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
#if defined(USE_SSE2) || defined(USE_AVX2)
#include <immintrin.h>
#endif
#include <math.h>
#include <assert.h>
#include <pthread.h>
#include "htslib/sam.h"
#include "../sigproc-flags.h"
#include "tailseq-dedup-approx.h"

/* Set a temporary limit of clusters in buffer. This will be
 * removed again once I implement a scalable routine that works for
 * larger clusters. */
#define XXX_TAGCLUSTER_LIMIT            4096

#define DEFAULT_TAGALN_BUFFER_LEN       8192
#define DEFAULT_TAGCLUSTER_BUFFER_LEN   8192

static inline int16_t
calculate_tag_prority(int flags)
{
#define SET(value, flagbit) (((value) & (flagbit)) != 0)
#define NOTSET(value, flagbit) (((value) & (flagbit)) == 0)
    return
        NOTSET(flags, PAFLAG_BARCODE_HAS_MISMATCHES) * 1024 +
        NOTSET(flags, PAFLAG_DELIMITER_NOT_FOUND) * 512 +
        NOTSET(flags, PAFLAG_DELIMITER_HAS_MISMATCH) * 256 +
        SET(flags, PAFLAG_MEASURED_FROM_FLUORESCENCE) * 128 +
        NOTSET(flags, PAFLAG_DARKCYCLE_OVER_THRESHOLD) * 64 +
        NOTSET(flags, PAFLAG_BALANCER_BIASED) * 32 +
        NOTSET(flags, PAFLAG_BALANCER_SIGNAL_BAD) * 16 +
        NOTSET(flags, PAFLAG_DARKCYCLE_EXISTS) * 8 +
        NOTSET(flags, PAFLAG_DELIMITER_IS_SHIFTED) * 4 +
        NOTSET(flags, PAFLAG_NO_POLYA_DETECTED) * 2 +
        NOTSET(flags, PAFLAG_MEASURED_USING_NAIVE_RULER);
#undef SET
#undef NOTSET
}

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
count_trimers(const char *seq, union trimer_composition *counts)
{
    /* This roughly translate trimer incidences to array 6-bit integers.
     * All non-ACGT characters are translated to 4, which may produce
     * inaccurate code leading to false positives. As this function is
     * used only to reduce comparisons between two very distant sequences,
     * false positives can be tolerated in favor of more frequently
     * evaluated instructions here. */
    int trimercode;
    const char *p;

    memset(counts->count, 0, 64);

//    fprintf(stderr, "---- %s -> ", seq);
    trimercode = (DNABASE2NUM[(int)seq[0]] << 2) | (DNABASE2NUM[(int)seq[1]]);
    for (p = seq + 2; *p != '\0'; p++) {
        trimercode = 63 & ((trimercode << 2) | DNABASE2NUM[(int)*p]);
//        fprintf(stderr, "%d ", trimercode);
        counts->count[trimercode]++;
    }
//    fprintf(stderr, "\n");
}

static struct deduppool *
deduppool_init(ssize_t tagaln_capacity, ssize_t tagcluster_capacity,
               ssize_t cologroup_capacity, int umi_length)
{
    struct deduppool *pool;

    pool = malloc(sizeof(struct deduppool));
    if (pool == NULL)
        return NULL;

    pool->tagalns = calloc(tagaln_capacity, sizeof(struct tagaln));
    if (pool->tagalns == NULL) {
        free(pool);
        return NULL;
    }
    pool->tagaln_capacity = tagaln_capacity;
    pool->tagaln_tryalloc_next = 0;

    pool->tagclusters = calloc(tagcluster_capacity, sizeof(struct tagcluster));
    if (pool->tagclusters == NULL) {
        free(pool->tagalns);
        free(pool);
        return NULL;
    }
    pool->tagcluster_capacity = tagcluster_capacity;
    pool->tagcluster_tryalloc_next = 0;

    pool->tagcluster_left = pool->tagcluster_right = -1;
    pool->tagclusters_live = 0;
    pool->umi_length = umi_length;

    memset(pool->tagalns, 0, sizeof(struct tagaln) * tagaln_capacity);
    memset(pool->tagclusters, 0, sizeof(struct tagcluster) * tagcluster_capacity);

    return pool;
}

static void
deduppool_destroy(struct deduppool *pool)
{
    free(pool->tagclusters);
    free(pool->tagalns);
    free(pool);
}

static ssize_t
tagaln_new(struct deduppool *pool)
{
    ssize_t i;

    for (i = pool->tagaln_tryalloc_next; i < pool->tagaln_capacity; i++)
        if (!is_tagaln_occupied(&pool->tagalns[i])) {
            pool->tagaln_tryalloc_next = (pool->tagaln_tryalloc_next + 1) % pool->tagaln_capacity;
            return i;
        }

    for (i = 0; i < pool->tagaln_tryalloc_next; i++)
        if (!is_tagaln_occupied(&pool->tagalns[i])) {
            pool->tagaln_tryalloc_next = (pool->tagaln_tryalloc_next + 1) % pool->tagaln_capacity;
            return i;
        }

    {
        /* Resize the pool to create spaces. */
        ssize_t newcapacity=pool->tagaln_capacity * 2;
        struct tagaln *newptr;

        newptr = realloc(pool->tagalns, sizeof(struct tagaln) * newcapacity);
        if (newptr == NULL)
            return -1;

        /* Empty newly allocated elements */
        memset(&newptr[pool->tagaln_capacity], 0,
               sizeof(struct tagaln) * pool->tagaln_capacity);

        pool->tagalns = newptr;
        pool->tagaln_tryalloc_next = pool->tagaln_capacity;
        pool->tagaln_capacity = newcapacity;

        return pool->tagaln_tryalloc_next++;
    }
}

static inline ssize_t
tagaln_create(struct deduppool *pool, const char *qname, uint32_t flags,
              uint32_t ndups, int16_t polya_len_1, int16_t polya_len_2)
{
    ssize_t idx;
    struct tagaln *aln;

    idx = tagaln_new(pool);
    if (idx < 0)
        return -1;

    aln = &pool->tagalns[idx];
    strcpy(aln->qname, qname);
    aln->flags = flags;
    aln->ndups = ndups;
    aln->polya_len_1 = polya_len_1;
    aln->polya_len_2 = polya_len_2;
    aln->priority = calculate_tag_prority(flags);
    aln->next = -1;

    return idx;
}

ssize_t
tagcluster_new(struct deduppool *pool)
{
    ssize_t i;

    for (i = pool->tagcluster_tryalloc_next; i < pool->tagcluster_capacity; i++)
        if (!is_tagcluster_occupied(&pool->tagclusters[i])) {
            pool->tagcluster_tryalloc_next = (pool->tagcluster_tryalloc_next + 1)
                                             % pool->tagcluster_capacity;
            return i;
        }

    for (i = 0; i < pool->tagcluster_tryalloc_next; i++)
        if (!is_tagcluster_occupied(&pool->tagclusters[i])) {
            pool->tagcluster_tryalloc_next = (pool->tagcluster_tryalloc_next + 1)
                                             % pool->tagcluster_capacity;
            return i;
        }

    {
        /* Resize the pool to create spaces. */
        ssize_t newcapacity=pool->tagcluster_capacity * 2;
        struct tagcluster *newptr;

        newptr = realloc(pool->tagclusters, sizeof(struct tagcluster) * newcapacity);
        if (newptr == NULL)
            return -1;

        /* Empty newly allocated elements */
        memset(&newptr[pool->tagcluster_capacity], 0,
               sizeof(struct tagcluster) * pool->tagcluster_capacity);

        pool->tagclusters = newptr;
        pool->tagcluster_tryalloc_next = pool->tagcluster_capacity;
        pool->tagcluster_capacity = newcapacity;

        return pool->tagcluster_tryalloc_next++;
    }
}

static inline ssize_t
tagcluster_create(struct deduppool *pool, int32_t tid, int32_t pos,
                  const char *umiseq, uint32_t umi_rep_ndups,
                  ssize_t tagaln)
{
    ssize_t idx;
    struct tagcluster *clstr;

    idx = tagcluster_new(pool);
    if (idx < 0)
        return -1;

    clstr = &pool->tagclusters[idx];
    strcpy(clstr->umi_rep, umiseq);
    count_trimers(umiseq, &clstr->trimer_counts);
    clstr->tid = tid;
    clstr->pos = pos;
    clstr->tagaln_head = clstr->tagaln_tail = tagaln;
    clstr->prev = -1;
    clstr->next = -1;

    return idx;
}

void
tagcluster_free(struct deduppool *pool, struct tagcluster *clstr)
{
    while (clstr->tagaln_head >= 0) {
        struct tagaln *t=&pool->tagalns[clstr->tagaln_head];
        clstr->tagaln_head = t->next;
        tagaln_free(t);
    }

    clstr->umi_rep[0] = '\0';
}

static struct tagcluster *
tagcluster_pushright(struct deduppool *pool, ssize_t idx)
{
    if (pool->tagcluster_right < 0)
        pool->tagcluster_left = pool->tagcluster_right = idx;
    else {
        struct tagcluster *right=&pool->tagclusters[pool->tagcluster_right];
        assert(right->next == -1);
        pool->tagclusters[idx].prev = pool->tagcluster_right;
        pool->tagcluster_right = right->next = idx;
    }

    pool->tagclusters_live++;
    return &pool->tagclusters[idx];
}

static inline struct tagcluster *
tagcluster_peekleft(struct deduppool *pool)
{
    if (pool->tagcluster_left < 0)
        return NULL;

    assert(pool->tagcluster_right >= 0);
    return &pool->tagclusters[pool->tagcluster_left];
}

static inline struct tagcluster *
tagcluster_popleft_not_after_pos(struct deduppool *pool, int32_t tid, int32_t pos)
{
    struct tagcluster *clstr;

    if (pool->tagcluster_left < 0)
        return NULL;

    assert(pool->tagcluster_right >= 0);
    clstr = &pool->tagclusters[pool->tagcluster_left];
    if (clstr->tid == tid && clstr->pos >= pos)
        return NULL;

    pool->tagcluster_left = clstr->next;
    if (clstr->next < 0)
        pool->tagcluster_right = -1;
    clstr->prev = clstr->next = -1;

    pool->tagclusters_live--;

    return clstr;
}

static inline struct tagcluster *
tagcluster_popleft(struct deduppool *pool)
{
    struct tagcluster *clstr;

    if (pool->tagcluster_left < 0)
        return NULL;

    assert(pool->tagcluster_right >= 0);
    clstr = &pool->tagclusters[pool->tagcluster_left];
    pool->tagcluster_left = clstr->next;
    if (clstr->next < 0)
        pool->tagcluster_right = -1;
    clstr->prev = clstr->next = -1;

    pool->tagclusters_live--;

    return clstr;
}

struct tagcluster *
tagcluster_popmiddle(struct deduppool *pool, struct tagcluster *clstr)
{
    assert(clstr->prev >= 0);

    pool->tagclusters[clstr->prev].next = clstr->next;
    if (clstr->next >= 0)
        pool->tagclusters[clstr->next].prev = clstr->prev;
    else
        pool->tagcluster_right = clstr->prev;
    clstr->prev = clstr->next = -1;

    pool->tagclusters_live--;

    return clstr;
}

static uint64_t
diffcount_trimer_compositions(union trimer_composition *a, union trimer_composition *b)
{
#if defined(USE_AVX2)
    union {
        __m256i v256i[2];
        uint64_t v64i[8];
    } res;

    res.v256i[0] = _mm256_sad_epu8(a->avxreg[0], b->avxreg[0]);
    res.v256i[1] = _mm256_sad_epu8(a->avxreg[1], b->avxreg[1]);

    return res.v64i[0] + res.v64i[1] + res.v64i[2] + res.v64i[3] +
           res.v64i[4] + res.v64i[5] + res.v64i[6] + res.v64i[7];
#elif defined(USE_SSE2)
    union {
        __m128i v128i[4];
        uint64_t v64i[8];
    } res;

    res.v128i[0] = _mm_sad_epu8(a->sse2reg[0], b->sse2reg[0]);
    res.v128i[1] = _mm_sad_epu8(a->sse2reg[1], b->sse2reg[1]);
    res.v128i[2] = _mm_sad_epu8(a->sse2reg[2], b->sse2reg[2]);
    res.v128i[3] = _mm_sad_epu8(a->sse2reg[3], b->sse2reg[3]);

    return res.v64i[0] + res.v64i[1] + res.v64i[2] + res.v64i[3] +
           res.v64i[4] + res.v64i[5] + res.v64i[6] + res.v64i[7];
#else
    int res, i;

    for (i = 0, res = 0; i < 64; i++)
        res += abs(a->count[i] - b->count[i]);

    return res;
#endif
}

static inline int
MIN3(int a, int b, int c)
{
    if (a < b)
        if (a < c)
            return a;
        else
            return c;
    else
        if (b < c)
            return b;
        else
            return c;
}

static int
edit_distance(const char *s1, const char *s2, int length, int maximum_distance)
{
    int x, y, lastdiag, olddiag, left, right;
    int column[length+1];

    left = 1;
    right = maximum_distance + 1;

    for (y = left; y <= right; y++)
        column[y] = y;

    for (x = 1; x <= length; x++) {
        int nleft, nright;

        column[0] = x;
        lastdiag = x - 1;
        nleft = nright = -1;

        for (y = left; y <= right; y++) {
            olddiag = column[y];
            column[y] = MIN3(column[y] + 1, column[y-1] + 1, lastdiag + (s1[y-1] == s2[x-1] ? 0 : 1));
            lastdiag = olddiag;
            if (column[y] <= maximum_distance) {
                if (nleft < 0)
                    nleft = y;
                nright = y;
            }
        }

        if (nleft < 0)
            return maximum_distance + 1;

        left = nleft;
        right = (nright >= length ? length : nright + 1);
    }

    return column[length];
}

static int
write_tagcluster(struct deduppool *pool, struct tagcluster *clstr,
                 pthread_mutex_t *lock)
{
#define polya_len_ref       polya_len_2
#define mean_polya_len_ref  mean_polya_len_2
    struct tagaln *tag;
    ssize_t tag_ix, bestprio_tag_ix, bestmatching_ix;
    int64_t total_dups, bestprio_polya_sum_1, bestprio_polya_sum_2, bestprio_polya_pos_dups;
    float mean_polya_len_1, mean_polya_len_2, bestmatching_residual;
    int best_priority, bestprio_tags;

//    fprintf(stderr, "----- FLUSHING group %d %s --------------\n", (int)group_ix, grp->umi_rep);
//    for (tag_ix = grp->tagaln; tag_ix >= 0; tag_ix = tag->next) {
//        tag = &pool->tagalns[tag_ix];
//        tag->priority = calculate_tag_prority(tag->flags);
//        fprintf(stderr, " TAG[%d]:%d = Za:%d ZD:%d\n", (int)tag_ix, (int)tag->priority,
//                    tag->polya_len, tag->ndups);
//    }
//    fprintf(stderr, "-------------------\n");

    best_priority = -1;
    bestprio_tag_ix = -1;
    bestprio_polya_sum_1 = bestprio_polya_sum_2 = bestprio_polya_pos_dups = 0;
    total_dups = 0;
    bestprio_tags = 0;

    for (tag_ix = clstr->tagaln_head; tag_ix >= 0; tag_ix = tag->next) {
        tag = &pool->tagalns[tag_ix];
        total_dups += tag->ndups;
        tag->priority = calculate_tag_prority(tag->flags);
        if (tag->priority < best_priority)
            continue;

        if (tag->priority > best_priority) {
            best_priority = tag->priority;
            bestprio_tag_ix = tag_ix;
            bestprio_tags = 1;

            if (tag->polya_len_ref > 0) {
                bestprio_polya_sum_1 = tag->polya_len_1 * tag->ndups;
                bestprio_polya_sum_2 = tag->polya_len_2 * tag->ndups;
                bestprio_polya_pos_dups = tag->ndups;
            }
            else
                bestprio_polya_sum_1 = bestprio_polya_sum_2 = bestprio_polya_pos_dups = 0;
        }
        else {
            bestprio_tags++;
            if (tag->polya_len_ref > 0) {
                bestprio_polya_sum_1 += tag->polya_len_1 * tag->ndups;
                bestprio_polya_sum_2 += tag->polya_len_2 * tag->ndups;
                bestprio_polya_pos_dups += tag->ndups;
            }
        }
    }

    if (bestprio_polya_pos_dups <= 0 || bestprio_tags <= 1) {
        /* If a single tag is the best, show it as a final representation. */
        pthread_mutex_lock(lock);
        printf("%d\t%d\t%ld\t%s\n", pool->tagalns[bestprio_tag_ix].polya_len_1,
               pool->tagalns[bestprio_tag_ix].polya_len_2,
               total_dups, pool->tagalns[bestprio_tag_ix].qname);
        pthread_mutex_unlock(lock);

        return 0;
    }

    /* Choose the most similar tag to the mean poly(A) length having the highest
     * priority. */
    mean_polya_len_1 = (float)bestprio_polya_sum_1 / (float)bestprio_polya_pos_dups;
    mean_polya_len_2 = (float)bestprio_polya_sum_2 / (float)bestprio_polya_pos_dups;
//    fprintf(stderr, " -- poly(A) mean: %.2f nt.\n", mean_polya_len);

    bestmatching_ix = bestprio_tag_ix;
    bestmatching_residual = fabsf(mean_polya_len_ref -
                                  pool->tagalns[bestprio_tag_ix].polya_len_ref);

//    fprintf(stderr, " -- bestprio(%d) residual = %.2f nt\n", (int)bestmatching_ix,
//            bestmatching_residual);

    for (tag_ix = clstr->tagaln_head; tag_ix >= 0; tag_ix = tag->next) {
        float residual;

        tag = &pool->tagalns[tag_ix];
        if (tag->priority != best_priority)
            continue;

        residual = fabsf(mean_polya_len_ref - tag->polya_len_ref);
        if (residual < bestmatching_residual) {
            bestmatching_ix = tag_ix;
            bestmatching_residual = residual;
        }
    }

    pthread_mutex_lock(lock);
    printf("%d\t%d\t%ld\t%s\n",
           mean_polya_len_1 >= 0.f ? (int)roundf(mean_polya_len_1) : -1,
           mean_polya_len_2 >= 0.f ? (int)roundf(mean_polya_len_2) : -1,
           total_dups, pool->tagalns[bestmatching_ix].qname);
    pthread_mutex_unlock(lock);

//    fprintf(stderr, " -- (final) bestprio(%d) residual = %.2f nt\n", (int)bestmatching_ix,
//            bestmatching_residual);

    return 0;
}

static int
merge_similar_umi_clusters(struct deduppool *tpool, struct tagcluster *query,
                           int editdist_threshold, int trimercomp_threshold)
{
    struct tagcluster *target;
    ssize_t target_ix;

    assert(query->next >= 0);
    target = query;

    while (target->next >= 0) {
        int dist;

        target_ix = target->next;
        target = &tpool->tagclusters[target_ix];

        if ((int)diffcount_trimer_compositions(&query->trimer_counts,
                &target->trimer_counts) >= trimercomp_threshold)
            dist = editdist_threshold;
        else
            dist = edit_distance(query->umi_rep, target->umi_rep, tpool->umi_length,
                                 editdist_threshold - 1);

//        fprintf(stderr, "[%d] DISTANCE %s-%s = %d\n", (int)target_ix,
//                query->umi_rep, target->umi_rep, dist);

        if (dist < editdist_threshold) {
            /* Merge the target into the query cluster. */
            if (query->umi_rep_ndups < target->umi_rep_ndups) {
                /* Change the cluster rep UMI sequence when target
                 * is more abundant. */
                strcpy(query->umi_rep, target->umi_rep);
                query->umi_rep_ndups = target->umi_rep_ndups;
                memcpy(&query->trimer_counts.count,
                       &target->trimer_counts.count, sizeof(query->trimer_counts));
            }

            if (query->tagaln_tail < 0) {
                query->tagaln_head = target->tagaln_head;
                query->tagaln_tail = target->tagaln_tail;
            }
            else if (target->tagaln_tail < 0)
                /* pass */;
            else {
                assert(tpool->tagalns[query->tagaln_tail].next == -1);
                tpool->tagalns[query->tagaln_tail].next = target->tagaln_head;
                query->tagaln_tail = target->tagaln_tail;
            }

            target->tagaln_head = target->tagaln_tail = -1;
            (void)tagcluster_popmiddle(tpool, target);
            tagcluster_free(tpool, target);
            target = query;
        }
    }

    return 0;
}

static int
perform_tag_clustering(struct deduppool *tpool, int editdist_threshold,
                       int trimercomp_threshold)
{
    struct tagcluster *clstr;

    if (tpool->tagclusters_live <= 1)
        return 0;

//    fprintf(stderr, "==> Tag clustering! CLUSTERS %d-%d // %d->",
//            (int)tpool->tagcluster_left,
//            (int)tpool->tagcluster_right,
//            (int)tpool->tagclusters_live);
    for (clstr = &tpool->tagclusters[tpool->tagcluster_left];
            clstr->next >= 0; clstr = &tpool->tagclusters[clstr->next]) {
        if (merge_similar_umi_clusters(tpool, clstr, editdist_threshold,
                                       trimercomp_threshold) < 0)
            return -1;

        /* Next pointer can be modified during cluster merging. */
        if (clstr->next < 0)
            break;
    }

    //fprintf(stderr, "%d\n", (int)tpool->tagclusters_live);
    return 0;
}

static int
deduplicate_tailseq_bam(struct taskpool *tasks, samFile *samf, hts_idx_t *bamidx,
                        struct deduppool *tpool, int tid)
{
    bam1_t *bamentry;
    hts_itr_t *bamiter;
    int coorddist_tolerance, editdist_threshold, trimercomp_threshold;
    int ret;

    coorddist_tolerance = tasks->coorddist_tolerance;
    editdist_threshold = tasks->editdist_tolerance + 1;
    trimercomp_threshold = tasks->editdist_tolerance * 6 + 1;

    bamentry = bam_init1();
    if (bamentry == NULL) {
        perror("bam_init1");
        return -1;
    }

    bamiter = bam_itr_queryi(bamidx, tid, 0, tasks->target_lengths[tid]);
    if (bamiter == NULL) {
        perror("bam_itr_queryi");
        return -1;
    }

    while ((ret = bam_itr_next(samf, bamiter, bamentry)) >= 0) {
        struct tagaln *newtag;
        uint8_t *umiseq_aux, *flags_aux, *palen1_aux, *palen2_aux, *ndups_aux;
        ssize_t tag_ix, cluster_ix;

        if (bamentry->core.flag & BAM_FUNMAP)
            continue;

        umiseq_aux = bam_aux_get(bamentry, "ZM");
        flags_aux = bam_aux_get(bamentry, "ZF");
        palen1_aux = bam_aux_get(bamentry, "ZA");
        palen2_aux = bam_aux_get(bamentry, "Za");
        ndups_aux = bam_aux_get(bamentry, "ZD");
        if (umiseq_aux == NULL || flags_aux == NULL || palen1_aux == NULL ||
                palen2_aux == NULL || ndups_aux == NULL)
            continue; /* Keep non-TAIL-seq entries untouched */

        if (bamentry->core.l_qname > QNAME_LEN_MAX) {
            fprintf(stderr, "ERROR: A query name too long \"%s\".\n",
                    (const char *)bamentry->data);
            return -1;
        }

        /* Flush out colocated groups became too far from the cursor. */
        if (!is_deduppool_empty(tpool)) {
            struct tagcluster *clstr;
            int32_t validfrom=bamentry->core.pos - coorddist_tolerance;
            int XXX_force_flushing=(tpool->tagclusters_live >= XXX_TAGCLUSTER_LIMIT);

            clstr = tagcluster_peekleft(tpool);
            /*if (clstr->tid != bamentry->core.tid || clstr->pos < validfrom) { XXX */
            if (XXX_force_flushing ||
                    clstr->tid != bamentry->core.tid || clstr->pos < validfrom) {
                if (perform_tag_clustering(tpool, editdist_threshold,
                                           trimercomp_threshold) < 0)
                    return -1;

                while ((XXX_force_flushing && (clstr = tagcluster_popleft(tpool)) != NULL) ||
                       (clstr = tagcluster_popleft_not_after_pos(tpool, bamentry->core.tid,
                                                                 validfrom)) != NULL) {
                    if (write_tagcluster(tpool, clstr, &tasks->writelock) < 0)
                        return -1;
                    tagcluster_free(tpool, clstr);
                }
            }
        }

        /* Create a tag and a cluster entries. */
        tag_ix = tagaln_create(tpool, (const char *)bamentry->data,
                               bam_aux2i(flags_aux), bam_aux2i(ndups_aux),
                               bam_aux2i(palen1_aux), bam_aux2i(palen2_aux));
        if (tag_ix < 0) {
            perror("tagaln_create");
            return -1;
        }
        newtag = &tpool->tagalns[tag_ix];

        cluster_ix = tagcluster_create(tpool, bamentry->core.tid, bamentry->core.pos,
                                       bam_aux2Z(umiseq_aux), newtag->ndups, tag_ix);
        if (cluster_ix < 0) {
            perror("tagcluster_create");
            return -1;
        }
        (void)tagcluster_pushright(tpool, cluster_ix);

//        fprintf(stderr, "\ntid=%d pos=%d qname=%s ZM=%s Za=%d ZF=%d ZD=%d\n",
//                (int)bamentry->core.tid, (int)bamentry->core.pos,
//                newtag->qname, newcluster->umi_rep, (int)newtag->polya_len_1,
//                (int)newtag->flags, (int)newtag->ndups);
    }

    if (!is_deduppool_empty(tpool)) {
        struct tagcluster *clstr;

        if (perform_tag_clustering(tpool, editdist_threshold,
                                   trimercomp_threshold) < 0)
            return -1;

        while ((clstr = tagcluster_popleft(tpool)) != NULL) {
            if (write_tagcluster(tpool, clstr, &tasks->writelock) < 0)
                return -1;
            tagcluster_free(tpool, clstr);
        }
    }

    bam_itr_destroy(bamiter);
    bam_destroy1(bamentry);

    return 0;
}

static int
get_bam_umi_length(samFile *samf, bam_hdr_t *header)
{
    bam1_t *aln=bam_init1();
    uint8_t *umiseq_aux;
    int r;

    if (aln == NULL)
        return -1;

    if (sam_read1(samf, header, aln) < 0) {
        bam_destroy1(aln);
        return -1;
    }

    umiseq_aux = bam_aux_get(aln, "ZM");
    if (umiseq_aux == NULL) {
        bam_destroy1(aln);
        return -1;
    }

    r = strlen(bam_aux2Z(umiseq_aux));

    bam_destroy1(aln);

    return r;
}

static void
process_deduplications(struct taskpool *tasks)
{
    samFile *samf;
    bam_hdr_t *header;
    hts_idx_t *bamidx;
    struct deduppool *tpool;

#define set_error_and_exit(code)    do { \
            pthread_mutex_lock(&tasks->poollock); \
            tasks->ret_code = (code); \
            pthread_mutex_unlock(&tasks->poollock); \
            pthread_exit((void *)0); \
        } while (0)

    samf = sam_open(tasks->bam_filename, "rb");
    if (samf == NULL) {
        fprintf(stderr, "ERROR: Can't open %s.\n", tasks->bam_filename);
        set_error_and_exit(1);
    }

    bamidx = sam_index_load(samf, tasks->bam_filename);
    if (bamidx == NULL) {
        fprintf(stderr, "ERROR: Failed to load the BAM index\n");
        set_error_and_exit(2);
    }

    header = sam_hdr_read(samf);
    if (header == NULL) {
        fprintf(stderr, "ERROR: Failed to read header from %s.\n", tasks->bam_filename);
        set_error_and_exit(3);
    }

    tpool = deduppool_init(DEFAULT_TAGALN_BUFFER_LEN, DEFAULT_TAGCLUSTER_BUFFER_LEN,
                           tasks->coorddist_tolerance + 2, tasks->umi_length);
    if (tpool == NULL) {
        perror("deduppool_init");
        set_error_and_exit(5);
    }

    for (;;) {
        int my_tid;

        pthread_mutex_lock(&tasks->poollock);
            if (tasks->tid_next > tasks->tid_max) {
                pthread_mutex_unlock(&tasks->poollock);
                break;
            }

            my_tid = tasks->tid_next++;
        pthread_mutex_unlock(&tasks->poollock);

        if (deduplicate_tailseq_bam(tasks, samf, bamidx, tpool, my_tid) < 0) {
            perror("deduplicate_tailseq_bam");
            set_error_and_exit(6);
        }
    }

    deduppool_destroy(tpool);

    hts_idx_destroy(bamidx);
    bam_hdr_destroy(header);

    sam_close(samf);

    pthread_exit((void *)0);
}

static int
load_sam_targets_count(struct taskpool *pool)
{
    samFile *samf;
    bam_hdr_t *header;
    hts_idx_t *bamidx;

    samf = sam_open(pool->bam_filename, "rb");
    if (samf == NULL) {
        fprintf(stderr, "ERROR: Can't open %s.\n", pool->bam_filename);
        return -1;
    }

    /* Try loading the index to fail early before the worker threads detect. */
    bamidx = sam_index_load(samf, pool->bam_filename);
    if (bamidx == NULL) {
        fprintf(stderr, "ERROR: Failed to load the BAM index\n");
        return -1;
    }

    header = sam_hdr_read(samf);
    if (header == NULL) {
        fprintf(stderr, "ERROR: Failed to read header from %s.\n", pool->bam_filename);
        return -1;
    }

    pool->umi_length = get_bam_umi_length(samf, header);
    if (pool->umi_length < 0) {
        fprintf(stderr, "ERROR: Can't verify the length of UMI from from %s.\n",
                pool->bam_filename);
        return -1;
    }

    pool->tid_max = header->n_targets - 1;
    pool->target_lengths = calloc(header->n_targets,
                                  sizeof(*header->target_len));
    memcpy(pool->target_lengths, header->target_len,
            sizeof(*header->target_len) * header->n_targets);

    hts_idx_destroy(bamidx);
    bam_hdr_destroy(header);

    sam_close(samf);

    return 0;
}

int
main(int argc, char *argv[])
{
    int nthreads;
    struct taskpool tasks;

    tasks.tid_next = 0;
    tasks.ret_code = 0;

    tasks.bam_filename = argv[1];
    tasks.coorddist_tolerance = atoi(argv[2]);
    tasks.editdist_tolerance = atoi(argv[3]);
    nthreads = atoi(argv[4]);

    pthread_mutex_init(&tasks.poollock, NULL);
    pthread_mutex_init(&tasks.writelock, NULL);

    if (load_sam_targets_count(&tasks) < 0)
        return 101;

    {
        pthread_t threads[nthreads];
        int i;

        for (i = 0; i < nthreads; i++)
            pthread_create(&threads[i], NULL, (void *)process_deduplications, (void *)&tasks);

        for (i = 0; i < nthreads; i++)
            pthread_join(threads[i], NULL);
    }

    pthread_mutex_destroy(&tasks.writelock);
    pthread_mutex_destroy(&tasks.poollock);

    return tasks.ret_code;
}
