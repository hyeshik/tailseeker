/*
 * tailseq-dedup-approx.h
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

#ifndef _TAILSEQ_DEDUP_APPROX_H_
#define _TAILSEQ_DEDUP_APPROX_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define QNAME_LEN_MAX                   32
#define UMI_LEN_MAX                     32

union trimer_composition {
    unsigned char count[64];
#ifdef USE_SSE2
    __m128i sse2reg[4];
#endif 
#ifdef USE_AVX2
    __m256i avxreg[2];
#endif
};

struct taskpool {
    int tid_next;
    int tid_max;
    uint32_t *target_lengths; /* borrowed pointer from sam header */

    int ret_code;

    const char *bam_filename;
    int coorddist_tolerance;
    int editdist_tolerance;
    int cdhit_bypass;

    pthread_mutex_t poollock;
    pthread_mutex_t writelock;
    pthread_mutex_t cdhitlock;
};

struct tagaln;
struct tagcluster;

#define is_tagaln_occupied(t)       ((t)->qname[0] != '\0')
#define tagaln_free(t)              do { (t)->qname[0] = '\0'; } while (0)

struct tagaln {
    char qname[QNAME_LEN_MAX];
    uint32_t flags;
    uint32_t ndups;
    int16_t polya_len;
    int16_t priority;
    ssize_t next; /* in a tagcluster */
    ssize_t __pad;
};

#define is_tagcluster_occupied(t)     ((t)->umi_rep[0] != '\0')

struct tagcluster {
    int32_t tid;
    int32_t pos;
    char umi_rep[UMI_LEN_MAX];
    uint32_t umi_rep_ndups;
    union trimer_composition trimer_counts;
    ssize_t tagaln_head;
    ssize_t tagaln_tail;
    ssize_t prev;
    ssize_t next;
};

#define is_deduppool_empty(pool)   ((pool)->tagcluster_left == -1)

struct deduppool {
    struct tagaln *tagalns;
    ssize_t tagaln_capacity;
    ssize_t tagaln_tryalloc_next;

    struct tagcluster *tagclusters;
    ssize_t tagcluster_capacity;
    ssize_t tagcluster_tryalloc_next;

    ssize_t tagcluster_left;
    ssize_t tagcluster_right;
    ssize_t tagclusters_live;
};

/* cdhit-proxy.cc */
extern int init_cdhit_internals(int nthreads, double similarity_threshold);
extern int cdhit_cluster_minitags(struct deduppool *pool);

/* tailseq-dedup-approx.c */
extern ssize_t tagcluster_new(struct deduppool *pool);
extern void tagcluster_free(struct deduppool *pool, struct tagcluster *clstr);
extern struct tagcluster *tagcluster_popmiddle(struct deduppool *pool,
                                               struct tagcluster *clstr);


#ifdef __cplusplus
};
#endif

#endif
