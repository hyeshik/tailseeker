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
#include "tailseq-sigproc.h"


static int
spawn_writers(const char *writercmd, struct SampleInfo *barcodes)
{
    for (; barcodes != NULL; barcodes = barcodes->next) {
        char cmd[BUFSIZ];

        snprintf(cmd, BUFSIZ-1, writercmd, barcodes->name);
        barcodes->stream = popen(cmd, "w");
        if (barcodes->stream == NULL) {
            perror("spawn_writers");
            fprintf(stderr, "Failed to run a writer subprocess: %s\n", cmd);
            return -1;
        }
    }

    return 0;
}


static void
terminate_writers(struct SampleInfo *barcodes)
{
    for (; barcodes != NULL; barcodes = barcodes->next)
        if (barcodes->stream != NULL) {
            pclose(barcodes->stream);
            barcodes->stream = NULL;
        }
}


static int
initialize_cif_bcl_buffers(int ncycles, int blocksize, struct CIFData ***intensities,
                           struct BCLData ***basecalls)
{
    struct CIFData **cifdata;
    struct BCLData **bcldata;
    int i;

    cifdata = malloc(sizeof(struct CIFData *) * ncycles);
    if (cifdata == NULL) {
        perror("initialize_cif_bcl_buffers");
        return -1;
    }

    bcldata = malloc(sizeof(struct BCLData *) * ncycles);
    if (bcldata == NULL) {
        free(cifdata);
        perror("initialize_cif_bcl_buffers");
        return -1;
    }

    memset(cifdata, 0, sizeof(struct CIFData *) * ncycles);
    memset(bcldata, 0, sizeof(struct BCLData *) * ncycles);

    for (i = 0; i < ncycles; i++) {
        cifdata[i] = new_cif_data(blocksize);
        if (cifdata[i] == NULL)
            goto onError;

        bcldata[i] = new_bcl_data(blocksize);
        if (bcldata[i] == NULL)
            goto onError;
    }

    *intensities = cifdata;
    *basecalls = bcldata;

    return 0;

  onError:
    perror("initialize_cif_bcl_buffers");

    for (i = 0; i < ncycles; i++) {
        if (cifdata[i] != NULL)
            free_cif_data(cifdata[i]);
        if (bcldata[i] != NULL)
            free_bcl_data(bcldata[i]);
    }

    free(bcldata);
    free(cifdata);

    return -1;
}


static int
load_intensities_and_basecalls(struct CIFReader **cifreader, struct BCLReader **bclreader,
                               struct AlternativeCallInfo *altcalls, int ncycles,
                               int blocksize, struct CIFData **intensities,
                               struct BCLData **basecalls)
{
    int cycleno;

    for (; altcalls != NULL; altcalls = altcalls->next)
        if (load_alternative_calls(altcalls->reader, basecalls + altcalls->first_cycle,
                                   blocksize) == -1)
            return -1;

    for (cycleno = 0; cycleno < ncycles; cycleno++) {
        if (load_cif_data(cifreader[cycleno], intensities[cycleno], blocksize) == -1)
            return -1;

        if (bclreader[cycleno] == BCLREADER_OVERRIDDEN)
            /* do nothing */;
        else if (load_bcl_data(bclreader[cycleno], basecalls[cycleno], blocksize) == -1)
            return -1;
    }

    return 0;
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
    const char *writercmd="/usr/local/bin/bgzip -c > scratch/demux-sqi/%s_a1101.sqi.gz";
    /* XXX change this writer command */

    blocksize = cfg->read_buffer_entry_count;

    snprintf(msgprefix, BUFSIZ, "[%s%d] ", cfg->laneid, cfg->tile);

    printf("%sOpening writer subprocesses\n", msgprefix);
    if (spawn_writers(writercmd, cfg->samples) == -1)
        return -1;

    printf("%sOpening input sources\n", msgprefix);

    cifreader = NULL;
    bclreader = NULL;
    intensities = NULL;
    basecalls = NULL;

    if (open_alternative_calls_bundle(msgprefix, cfg->altcalls) == -1)
        return -1;

    cifreader = open_cif_readers(msgprefix, cfg->datadir, cfg->lane, cfg->tile,
                                 cfg->total_cycles);
    if (cifreader == NULL)
        goto onError;

    bclreader = open_bcl_readers(msgprefix, cfg->datadir, cfg->lane, cfg->tile,
                                 cfg->total_cycles, cfg->altcalls);
    if (bclreader == NULL)
        goto onError;

    if (initialize_cif_bcl_buffers(cfg->total_cycles, blocksize,
                                   &intensities, &basecalls) == -1)
        return -1;

    clusters_to_go = nclusters = cifreader[0]->nclusters;
    totalblocks = nclusters / blocksize + ((nclusters % blocksize > 0) ? 1 : 0);
    printf("%sProcessing %u clusters.\n", msgprefix, nclusters);

    for (blockno = 0; clusters_to_go > 0; blockno++) {
        clusters_to_read = (clusters_to_go >= blocksize) ? blocksize : clusters_to_go;

        snprintf(msgprefix, BUFSIZ, "[%s%d#%d/%d] ", cfg->laneid, cfg->tile, blockno + 1,
                                                     totalblocks);
        printf("%sLoading CIF and BCL files\n", msgprefix);

        if (load_intensities_and_basecalls(cifreader, bclreader, cfg->altcalls,
                                           cfg->total_cycles,
                                           clusters_to_read, intensities, basecalls) == -1)
            goto onError;

        printf("%sDemultiplexing and writing\n", msgprefix);
        if (process_spots(cfg, nclusters - clusters_to_go,
                          intensities, basecalls) < 0)
            goto onError;

        clusters_to_go -= clusters_to_read;
    }

    printf("%sClearing\n", msgprefix);
    for (cycleno = 0; cycleno < cfg->total_cycles; cycleno++) {
        free_cif_data(intensities[cycleno]);
        free_bcl_data(basecalls[cycleno]);
    }

    free(intensities);
    free(basecalls);

    close_bcl_readers(bclreader, cfg->total_cycles);
    close_cif_readers(cifreader, cfg->total_cycles);
    if (close_alternative_calls_bundle(cfg->altcalls, 1) < 0)
        goto onError;

    terminate_writers(cfg->samples);

    printf("[%s%d] Finished.\n", cfg->laneid, cfg->tile);

    return 0;

  onError:
    close_alternative_calls_bundle(cfg->altcalls, 0);

    if (cifreader != NULL)
        close_cif_readers(cifreader, cfg->total_cycles);

    if (bclreader != NULL)
        close_bcl_readers(bclreader, cfg->total_cycles);

    for (cycleno = 0; cycleno < cfg->total_cycles; cycleno++) {
        if (intensities != NULL && intensities[cycleno] != NULL)
            free_cif_data(intensities[cycleno]);

        if (basecalls != NULL && basecalls[cycleno] != NULL)
            free_bcl_data(basecalls[cycleno]);
    }

    free(intensities);
    free(basecalls);

    terminate_writers(cfg->samples);

    return -1;
}


static int
write_demultiplexing_statistics(const char *output, struct SampleInfo *barcodes)
{
    FILE *fp;

    fp = fopen(output, "w");
    if (fp == NULL)
        return -1;

    fprintf(fp, "name,index,max_index_mm,delim,delim_pos,max_delim_mm,"
                "cln_no_index_mm,cln_1_index_mm,cln_2+_index_mm,"
                "cln_fp_mm,cln_no_delim\n");

    for (; barcodes != NULL; barcodes = barcodes->next)
        fprintf(fp, "%s,%s,%d,%s,%d,%d,%d,%d,%d,%d,%d\n",
                barcodes->name, barcodes->index, barcodes->maximum_index_mismatches,
                barcodes->delimiter, barcodes->delimiter_pos,
                barcodes->maximum_delimiter_mismatches, barcodes->clusters_mm0,
                barcodes->clusters_mm1, barcodes->clusters_mm2plus,
                barcodes->clusters_fpmismatch, barcodes->clusters_nodelim);

    fclose(fp);

    return 0;
}


static void
usage(const char *prog)
{
    printf("\
tailseq-sigproc 3.0\
\n - collects intensities and base calls from Illumina sequencing for TAIL-seq\
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
