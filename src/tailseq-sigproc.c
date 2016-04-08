/*
 * tailseq-sigproc.c
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
#include "htslib/bgzf.h"
#include "tailseq-sigproc.h"


static int
open_writers(struct TailseekerConfig *cfg)
{
    struct SampleInfo *sample;

    for (sample = cfg->samples; sample != NULL; sample = sample->next) {
        char *filename, *fastqfile_format;

        fastqfile_format = replace_placeholder(cfg->fastq_output, "{name}",
                                               sample->name);

        filename = replace_placeholder(fastqfile_format, "{read}", "R5");
        sample->stream_fastq_5 = bgzf_open(filename, "w");
        if (sample->stream_fastq_5 == NULL) {
            perror("open_writers");
            fprintf(stderr, "Failed to write to %s\n", filename);
            free(filename);
            return -1;
        }
        free(filename);

        filename = replace_placeholder(fastqfile_format, "{read}", "R3");
        sample->stream_fastq_3 = bgzf_open(filename, "w");
        if (sample->stream_fastq_3 == NULL) {
            perror("open_writers");
            fprintf(stderr, "Failed to write to %s\n", filename);
            free(filename);
            return -1;
        }
        free(filename);

        free(fastqfile_format);

        if (cfg->taginfo_output != NULL) {
            filename = replace_placeholder(cfg->taginfo_output,
                                           "{name}", sample->name);
            sample->stream_taginfo = bgzf_open(filename, "w");
            if (sample->stream_taginfo == NULL) {
                perror("open_writers");
                fprintf(stderr, "Failed to write to %s\n", filename);
                free(filename);
                return -1;
            }
            free(filename);
        }

        if (cfg->signal_dump_output != NULL && sample->dump_processed_signals) {
            filename = replace_placeholder(cfg->signal_dump_output,
                                           "{name}", sample->name);
            sample->stream_signal_dump = bgzf_open(filename, "w");
            if (sample->stream_signal_dump == NULL) {
                perror("open_writers");
                fprintf(stderr, "Failed to write to %s\n", filename);
                free(filename);
                return -1;
            }
            free(filename);
        }
    }

    return 0;
}


static void
close_writers(struct SampleInfo *sample)
{
    for (; sample != NULL; sample = sample->next) {
        if (sample->stream_fastq_5 != NULL) {
            bgzf_close(sample->stream_fastq_5);
            sample->stream_fastq_5 = NULL;
        }

        if (sample->stream_fastq_3 != NULL) {
            bgzf_close(sample->stream_fastq_3);
            sample->stream_fastq_3 = NULL;
        }

        if (sample->stream_taginfo != NULL) {
            bgzf_close(sample->stream_taginfo);
            sample->stream_taginfo = NULL;
        }

        if (sample->stream_signal_dump != NULL) {
            bgzf_close(sample->stream_signal_dump);
            sample->stream_signal_dump = NULL;
        }
    }
}


static int
initialize_cif_bcl_buffers(struct TailseekerConfig *cfg,
                           struct CIFData ***intensities, struct BCLData ***basecalls)
{
    struct CIFData **cifdata;
    struct BCLData **bcldata;
    int i;

    cifdata = malloc(sizeof(struct CIFData *) * cfg->threep_length);
    if (cifdata == NULL) {
        perror("initialize_cif_bcl_buffers");
        return -1;
    }

    bcldata = malloc(sizeof(struct BCLData *) * cfg->total_cycles);
    if (bcldata == NULL) {
        free(cifdata);
        perror("initialize_cif_bcl_buffers");
        return -1;
    }

    memset(cifdata, 0, sizeof(struct CIFData *) * cfg->threep_length);
    memset(bcldata, 0, sizeof(struct BCLData *) * cfg->total_cycles);

    for (i = 0; i < cfg->threep_length; i++) {
        cifdata[i] = new_cif_data(cfg->read_buffer_entry_count);
        if (cifdata[i] == NULL)
            goto onError;
    }

    for (i = 0; i < cfg->total_cycles; i++) {
        bcldata[i] = new_bcl_data(cfg->read_buffer_entry_count);
        if (bcldata[i] == NULL)
            goto onError;
    }

    *intensities = cifdata;
    *basecalls = bcldata;

    return 0;

  onError:
    perror("initialize_cif_bcl_buffers");

    for (i = 0; i < cfg->threep_length; i++)
        if (cifdata[i] != NULL)
            free_cif_data(cifdata[i]);

    for (i = 0; i < cfg->total_cycles; i++)
        if (bcldata[i] != NULL)
            free_bcl_data(bcldata[i]);

    free(bcldata);
    free(cifdata);

    return -1;
}


static int
load_intensities_and_basecalls(struct TailseekerConfig *cfg,
                               struct CIFReader **cifreader, struct BCLReader **bclreader,
                               int blocksize,
                               struct CIFData **intensities, struct BCLData **basecalls)
{
    struct AlternativeCallInfo *altcalls;
    int cycleno;

    for (altcalls = cfg->altcalls; altcalls != NULL; altcalls = altcalls->next)
        if (load_alternative_calls(altcalls->reader, basecalls + altcalls->first_cycle,
                                   blocksize) == -1)
            return -1;

    for (cycleno = 0; cycleno < cfg->threep_length; cycleno++)
        if (load_cif_data(cifreader[cycleno], intensities[cycleno], blocksize) == -1)
            return -1;

    for (cycleno = 0; cycleno < cfg->total_cycles; cycleno++)
        if (bclreader[cycleno] == BCLREADER_OVERRIDDEN)
            /* do nothing */;
        else if (load_bcl_data(bclreader[cycleno], basecalls[cycleno], blocksize) == -1)
            return -1;

    return 0;
}


static struct ParallelJobPool *
prepare_split_jobs(struct SampleInfo *samples, uint32_t nclusters, int clusters_per_job)
{
    struct ParallelJobPool *pool;
    size_t poolmemsize;
    int njobs, i;

    njobs = (nclusters + clusters_per_job - 1) / clusters_per_job;
    poolmemsize = sizeof(struct ParallelJobPool) + sizeof(struct ParallelJob) * njobs;
    pool = malloc(poolmemsize);
    if (pool == NULL)
        return NULL;

    pool->job_next = pool->jobs_done = 0;
    pool->jobs_total = njobs;
    pthread_mutex_init(&pool->poollock, NULL);

    for (i = 0; i < njobs; i++) {
        pool->jobs[i].jobid = i;
        pool->jobs[i].start = i * clusters_per_job;
        pool->jobs[i].end = (i + 1) * clusters_per_job;
        if (pool->jobs[i].end > nclusters)
            pool->jobs[i].end = nclusters;
    }

    for (; samples != NULL; samples = samples->next) {
        samples->wsync_fastq_5.jobs_written = 0;
        samples->wsync_fastq_3.jobs_written = 0;
        samples->wsync_taginfo.jobs_written = 0;

        pthread_cond_init(&samples->wsync_fastq_5.wakeup, NULL);
        pthread_cond_init(&samples->wsync_fastq_3.wakeup, NULL);
        pthread_cond_init(&samples->wsync_taginfo.wakeup, NULL);

        pthread_mutex_init(&samples->wsync_fastq_5.lock, NULL);
        pthread_mutex_init(&samples->wsync_fastq_3.lock, NULL);
        pthread_mutex_init(&samples->wsync_taginfo.lock, NULL);
    }

    return pool;
}


static void
free_parallel_jobs(struct ParallelJobPool *pool, struct SampleInfo *samples)
{
    pthread_mutex_destroy(&pool->poollock);
    free(pool);

    for (; samples != NULL; samples = samples->next) {
        pthread_cond_destroy(&samples->wsync_fastq_5.wakeup);
        pthread_cond_destroy(&samples->wsync_fastq_3.wakeup);
        pthread_cond_destroy(&samples->wsync_taginfo.wakeup);

        pthread_mutex_destroy(&samples->wsync_fastq_5.lock);
        pthread_mutex_destroy(&samples->wsync_fastq_3.lock);
        pthread_mutex_destroy(&samples->wsync_taginfo.lock);
    }
}


static int
run_spot_processing(struct ParallelJobPool *pool)
{
    struct WriteBuffer *wbuf, *wbuf0;
    size_t memsize, wbufsize;
    char *buf, *buf0;
    int i, r=0;

    memsize = (pool->bufsize_fastq_5 + pool->bufsize_fastq_3 + pool->bufsize_taginfo) *
              pool->cfg->num_samples;
    buf = buf0 = malloc(memsize);
    if (buf == NULL)
        goto onError;

    wbufsize = sizeof(struct WriteBuffer) * pool->cfg->num_samples;
    wbuf = malloc(wbufsize);
    if (wbuf == NULL) {
        free(buf);
        goto onError;
    }

    wbuf0 = malloc(wbufsize);
    if (wbuf0 == NULL) {
        free(buf);
        free(wbuf);
        goto onError;
    }

    for (i = 0; i < pool->cfg->num_samples; i++) {
        wbuf0[i].buf_fastq_5 = buf;
        buf += pool->bufsize_fastq_5;
        wbuf0[i].buf_fastq_3 = buf;
        buf += pool->bufsize_fastq_3;
        wbuf0[i].buf_taginfo = buf;
        buf += pool->bufsize_taginfo;
    }

    while (1) {
        struct ParallelJob *job;

        { /* Select a job to run. */
            pthread_mutex_lock(&pool->poollock);

            if (pool->error_occurred > 0 || pool->job_next >= pool->jobs_total) {
                pthread_mutex_unlock(&pool->poollock);
                break;
            }
            job = &pool->jobs[pool->job_next++];

            pthread_mutex_unlock(&pool->poollock);
        }

        memcpy(wbuf, wbuf0, wbufsize);

        r = process_spots(pool->cfg, pool->firstclusterno, pool->intensities,
                          pool->basecalls, wbuf0, wbuf, job->jobid, job->start,
                          job->end);
        if (r < 0)
            break;

        pthread_mutex_lock(&pool->poollock);
        pool->jobs_done++;
        pthread_mutex_unlock(&pool->poollock);
    }

    free(wbuf0);
    free(wbuf);
    free(buf0);

    if (r >= 0) {
        pthread_exit((void *)0);
        return 0;
    }

onError:
    pthread_mutex_lock(&pool->poollock);
    pool->error_occurred++;
    pthread_mutex_unlock(&pool->poollock);

    pthread_exit((void *)1);
    return -1;
}


static int
distribute_processing(struct TailseekerConfig *cfg, struct CIFData **intensities,
                      struct BCLData **basecalls, uint32_t firstclusterno)
{
    uint32_t cycleno, clustersinblock;
    struct ParallelJobPool *pool;
    pthread_t threads[cfg->threads];
    int i, r;

    clustersinblock = intensities[0]->nclusters;

    for (cycleno = 0; cycleno < cfg->threep_length; cycleno++)
        if (clustersinblock != intensities[cycleno]->nclusters) {
            fprintf(stderr, "Inconsistent number of clusters in CIF cycle %d.\n",
                    cfg->threep_length + cycleno + 1);
            return -1;
        }
    for (cycleno = 0; cycleno < cfg->total_cycles; cycleno++)
        if (clustersinblock != basecalls[cycleno]->nclusters) {
            fprintf(stderr, "Inconsistent number of clusters in cycle %d.\n", cycleno + 1);
            return -1;
        }

    pool = prepare_split_jobs(cfg->samples, clustersinblock, NUM_CLUSTERS_PER_JOB);
    if (pool == NULL)
        return -1;

    pool->error_occurred = 0;
    pool->cfg = cfg;
    pool->intensities = intensities;
    pool->basecalls = basecalls;
    pool->firstclusterno = firstclusterno;
    pool->bufsize_fastq_5 = NUM_CLUSTERS_PER_JOB * cfg->max_bufsize_fastq_5;
    pool->bufsize_fastq_3 = NUM_CLUSTERS_PER_JOB * cfg->max_bufsize_fastq_3;
    pool->bufsize_taginfo = NUM_CLUSTERS_PER_JOB * cfg->max_bufsize_taginfo;

    for (i = 0; i < cfg->threads; i++)
        pthread_create(&threads[i], NULL, (void *)run_spot_processing, (void *)pool);

    for (i = 0; i < cfg->threads; i++)
        pthread_join(threads[i], NULL);

    r = (pool->error_occurred > 0) * -1;

    free_parallel_jobs(pool, cfg->samples);

    return r;
}


static int
process(struct TailseekerConfig *cfg)
{
    struct CIFReader **cifreader;
    struct BCLReader **bclreader;
    struct CIFData **intensities;
    struct BCLData **basecalls;
    uint32_t clusters_to_go, blockno, nclusters, totalblocks;
    int cycleno, clusters_to_read, blocksize;
    char msgprefix[BUFSIZ];

    blocksize = cfg->read_buffer_entry_count;
    initialize_control_aligner(&cfg->controlinfo);

    snprintf(msgprefix, BUFSIZ, "[%s%d] ", cfg->laneid, cfg->tile);

    if (open_writers(cfg) == -1)
        return -1;

    printf("%sOpening input sources\n", msgprefix);

    cifreader = NULL;
    bclreader = NULL;
    intensities = NULL;
    basecalls = NULL;

    if (open_alternative_calls_bundle(msgprefix, cfg->altcalls) == -1)
        return -1;

    cifreader = open_cif_readers(msgprefix, cfg->datadir, cfg->lane, cfg->tile,
                                 cfg->threep_start, cfg->threep_length);
    if (cifreader == NULL)
        goto onError;

    bclreader = open_bcl_readers(msgprefix, cfg->datadir, cfg->lane, cfg->tile,
                                 cfg->total_cycles, cfg->altcalls);
    if (bclreader == NULL)
        goto onError;

    if (initialize_cif_bcl_buffers(cfg, &intensities, &basecalls) == -1)
        return -1;

    clusters_to_go = nclusters = cifreader[0]->nclusters;
    totalblocks = nclusters / blocksize + ((nclusters % blocksize > 0) ? 1 : 0);
    printf("%sProcessing %u clusters.\n", msgprefix, nclusters);

    for (blockno = 0; clusters_to_go > 0; blockno++) {
        clusters_to_read = (clusters_to_go >= blocksize) ? blocksize : clusters_to_go;

        snprintf(msgprefix, BUFSIZ, "[%s%d#%d/%d] ", cfg->laneid, cfg->tile, blockno + 1,
                                                     totalblocks);
        printf("%sLoading CIF and BCL files\n", msgprefix);

        if (load_intensities_and_basecalls(cfg, cifreader, bclreader,
                                           clusters_to_read, intensities, basecalls) == -1)
            goto onError;

        printf("%sAnalyzing and writing out\n", msgprefix);

        if (distribute_processing(cfg, intensities, basecalls,
                                  nclusters - clusters_to_go) < 0)
            goto onError;

        clusters_to_go -= clusters_to_read;
    }

    printf("%sClearing\n", msgprefix);
    for (cycleno = 0; cycleno < cfg->threep_length; cycleno++)
        free_cif_data(intensities[cycleno]);
    for (cycleno = 0; cycleno < cfg->total_cycles; cycleno++)
        free_bcl_data(basecalls[cycleno]);

    free(intensities);
    free(basecalls);

    close_bcl_readers(bclreader, cfg->total_cycles);
    close_cif_readers(cifreader, cfg->threep_length);
    if (close_alternative_calls_bundle(cfg->altcalls, 1) < 0)
        goto onError;

    close_writers(cfg->samples);
    free_control_aligner(&cfg->controlinfo);

    printf("[%s%d] Finished.\n", cfg->laneid, cfg->tile);

    return 0;

  onError:
    close_alternative_calls_bundle(cfg->altcalls, 0);

    if (cifreader != NULL)
        close_cif_readers(cifreader, cfg->threep_length);

    if (bclreader != NULL)
        close_bcl_readers(bclreader, cfg->total_cycles);

    for (cycleno = 0; cycleno < cfg->threep_length; cycleno++)
        if (intensities != NULL && intensities[cycleno] != NULL)
            free_cif_data(intensities[cycleno]);

    for (cycleno = 0; cycleno < cfg->total_cycles; cycleno++)
        if (basecalls != NULL && basecalls[cycleno] != NULL)
            free_bcl_data(basecalls[cycleno]);

    free(intensities);
    free(basecalls);

    close_writers(cfg->samples);
    free_control_aligner(&cfg->controlinfo);

    return -1;
}


static int
write_demultiplexing_statistics(const char *output, struct SampleInfo *samples)
{
    FILE *fp;

    fp = fopen(output, "w");
    if (fp == NULL) {
        perror("write_demultiplexing_statistics");
        fprintf(stderr, "Failed to write the demultiplexing statistics: %s\n",
                output);
        return -1;
    }

    fprintf(fp, "name,index,max_index_mm,delim,delim_pos,max_delim_mm,"
                "cln_no_index_mm,cln_1_index_mm,cln_2+_index_mm,"
                "cln_fp_mm,cln_qc_fail,cln_no_delim\n");

    for (; samples != NULL; samples = samples->next)
        fprintf(fp, "%s,%s,%d,%s,%d,%d,%d,%d,%d,%d,%d,%d\n",
                samples->name, samples->index, samples->maximum_index_mismatches,
                samples->delimiter, samples->delimiter_pos,
                samples->maximum_delimiter_mismatches, samples->clusters_mm0,
                samples->clusters_mm1, samples->clusters_mm2plus,
                samples->clusters_fpmismatch, samples->clusters_qcfailed,
                samples->clusters_nodelim);

    fclose(fp);

    return 0;
}


static void
usage(const char *prog)
{
    printf("\
tailseq-sigproc 3.0\
\n - Processes Illumina .cif and .bcl files to produce poly(A) lengths and\
\n   additional 3' end modifications.\
\n\
\nUsage: %s {config.ini}\
\n\
\nMail bug reports and suggestions to Hyeshik Chang <hyeshik@snu.ac.kr>.\n\n", prog);
}


int
main(int argc, char *argv[])
{
    struct TailseekerConfig *cfg;
    int r;

    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }

    cfg = parse_config(argv[1]);
    if (cfg == NULL)
        return -1;

    if (load_color_matrix(cfg->rulerparams.colormatrix,
                          cfg->threep_colormatrix_filename) < 0)
        return -1;

    r = process(cfg);

    if (r == 0 && cfg->stats_output != NULL)
        r = write_demultiplexing_statistics(cfg->stats_output, cfg->samples);

    free_config(cfg);

    return r;
}
