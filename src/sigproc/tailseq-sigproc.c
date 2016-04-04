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


#define DEFAULT_BATCH_BLOCK_SIZE    262144
/* BCL - 2 bytes per cluster per cycle
 * CIF - 8 bytes per cluster per cycle
 * 3080 bytes per cluster for regular 51+251 set-up.
 * then, 262144 clusters takes 770 MiB up for storing data. */


static int
spawn_writers(const char *writercmd, struct BarcodeInfo *barcodes)
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
terminate_writers(struct BarcodeInfo *barcodes)
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
process(const char *datadir, const char *laneid, int lane, int tile, int ncycles,
        int threep_start, int threep_length, int index_start, int index_length,
        struct BarcodeInfo *barcodes, const char *writercmd,
        struct AlternativeCallInfo *altcalls, struct ControlFilterInfo *control_info,
        struct PolyAFinderParameters *finderparams,
        struct PolyARulerParameters *rulerparams, int blocksize, int keep_no_delimiter)
{
    struct CIFReader **cifreader;
    struct BCLReader **bclreader;
    struct CIFData **intensities;
    struct BCLData **basecalls;
    uint32_t clusters_to_go, blockno, nclusters, totalblocks;
    int cycleno, clusters_to_read;
    char msgprefix[BUFSIZ];

    snprintf(msgprefix, BUFSIZ, "[%s%d] ", laneid, tile);

    printf("%sOpening writer subprocesses\n", msgprefix);
    if (spawn_writers(writercmd, barcodes) == -1)
        return -1;

    printf("%sOpening input sources\n", msgprefix);

    cifreader = NULL;
    bclreader = NULL;
    intensities = NULL;
    basecalls = NULL;

    if (open_alternative_calls_bundle(msgprefix, altcalls) == -1)
        return -1;

    cifreader = open_cif_readers(msgprefix, datadir, lane, tile, ncycles);
    if (cifreader == NULL)
        goto onError;

    bclreader = open_bcl_readers(msgprefix, datadir, lane, tile, ncycles, altcalls);
    if (bclreader == NULL)
        goto onError;

    if (initialize_cif_bcl_buffers(ncycles, blocksize, &intensities, &basecalls) == -1)
        return -1;

    clusters_to_go = nclusters = cifreader[0]->nclusters;
    totalblocks = nclusters / blocksize + ((nclusters % blocksize > 0) ? 1 : 0);
    printf("%sProcessing %u clusters.\n", msgprefix, nclusters);

    for (blockno = 0; clusters_to_go > 0; blockno++) {
        clusters_to_read = (clusters_to_go >= blocksize) ? blocksize : clusters_to_go;

        snprintf(msgprefix, BUFSIZ, "[%s%d#%d/%d] ", laneid, tile, blockno + 1, totalblocks);
        printf("%sLoading CIF and BCL files\n", msgprefix);

        if (load_intensities_and_basecalls(cifreader, bclreader, altcalls, ncycles,
                                           clusters_to_read, intensities, basecalls) == -1)
            goto onError;

        printf("%sDemultiplexing and writing\n", msgprefix);
        if (process_spots(laneid, tile, ncycles, nclusters - clusters_to_go,
                          threep_start, threep_length, index_start, index_length,
                          barcodes, intensities, basecalls, control_info,
                          finderparams, rulerparams, keep_no_delimiter) == -1)
            goto onError;

        clusters_to_go -= clusters_to_read;
    }

    printf("%sClearing\n", msgprefix);
    for (cycleno = 0; cycleno < ncycles; cycleno++) {
        free_cif_data(intensities[cycleno]);
        free_bcl_data(basecalls[cycleno]);
    }

    free(intensities);
    free(basecalls);

    close_bcl_readers(bclreader, ncycles);
    close_cif_readers(cifreader, ncycles);
    if (close_alternative_calls_bundle(altcalls, 1) < 0)
        goto onError;

    terminate_writers(barcodes);

    printf("[%s%d] Finished.\n", laneid, tile);

    return 0;

  onError:
    close_alternative_calls_bundle(altcalls, 0);

    if (cifreader != NULL)
        close_cif_readers(cifreader, ncycles);

    if (bclreader != NULL)
        close_bcl_readers(bclreader, ncycles);

    for (cycleno = 0; cycleno < ncycles; cycleno++) {
        if (intensities != NULL && intensities[cycleno] != NULL)
            free_cif_data(intensities[cycleno]);

        if (basecalls != NULL && basecalls[cycleno] != NULL)
            free_bcl_data(basecalls[cycleno]);
    }

    free(intensities);
    free(basecalls);

    terminate_writers(barcodes);

    return -1;
}


static int
write_demultiplexing_statistics(const char *output, struct BarcodeInfo *barcodes)
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
\nUsage: %s [OPTION]...\
\n\
\nRequired parameters:\
\n  -d,  --data-dir=DIR               Illumina \"Intensities\" directory.\
\n  -r,  --run-id=ID                  short identifier of the run.\
\n  -l,  --lane=NUM                   lane number.\
\n  -t,  --tile=NUM                   tile number.\
\n  -n,  --ncycles=NUM                number of cycles.\
\n  -i,  --index-read=BEGIN,LENGTH    cycle number range of the index read.\
\n  -e,  --threep-read=BEGIN,LENGTH   cycle number range of the 3'-end read.\
\n  -w,  --writer-command=COMMAND     shell command to run to write .sqi\
\n                                    output (usually a stream compressor.)\
\n  -x,  --color-matrix=PATH          fluorescence crosstalk estimation matrix.\
\n\
\nOptional parameters:\
\n  -c,  --alternative-call=SPEC      FASTQ file to replace base calls.\
\n                                    specify as \"filename,first_cycle\".\
\n  -m,  --sample=SPEC                sample information as formatted in\
\n            \"name,index,max_index_miss,delimiter,delimiter_pos,max_delim_miss,fingerprint,max_fingerprint_miss,umi_pos,umi_length\".\
\n       --keep-no-delimiter          don't skip reads where a delimiter\
\n                                    is not found.\
\n  -f,  --filter-control=SPEC        sort out PhiX control reads by sequence\
\n                                    in \"name,first_cycle,length\".\
\n  -B,  --block-size=NUM             number of clusters in a processing unit\
\n                                    (default: 262144).\
\n  -o,  --demultiplex-stats=PATH     file path where write cluster count statistics\
\n                                    on demultiplexing (default: no output)\
\n  -s,  --sigproc-trigger=NUM        poly(A) length to handle over to fluorescence\
\n                                    signal processor (default: 10)\
\n       --minimum-polya-length=NUM   suppress poly(A) calling when shorter \
\n                                    than this. (default: 5)\
\n       --maximum-modifications=NUM  maximum length of non-poly(A) modifications \
\n                                    to poly(A). (default: 10)\
\n       --weight-T=NUM               poly(A) finder weight for T (default: 1) \
\n       --weight-ACG=NUM             poly(A) finder weight for A/C/G \
\n                                    (default: -10) \
\n       --weight-N=NUM               poly(A) finder weight for N (default: -5) \
\n\
\nAll cycle numbers or positions follow the 1-based inclusive system.\
\n\
\nMail bug reports and suggestions to Hyeshik Chang <hyeshik@snu.ac.kr>.\n\n", prog);
}


int
main(int argc, char *argv[])
{
    char *datadir, *runid, *writercmd;
    char *demultiplex_stats_filename;
    char *color_matrix_filename;
    int keep_no_delimiter_flag;
    int lane, tile, ncycles, blocksize;
    int index_start, index_length;
    int threep_start, threep_length;
    struct BarcodeInfo *barcodes;
    struct AlternativeCallInfo *altcalls;
    struct ControlFilterInfo controlinfo;
    struct PolyAFinderParameters finderparams;
    struct PolyARulerParameters rulerparams;
    int c, r;

    struct option long_options[] =
    {
        {"keep-no-delimiter",   no_argument,        &keep_no_delimiter_flag,    1},
        {"data-dir",            required_argument,  0,                          'd'},
        {"run-id",              required_argument,  0,                          'r'},
        {"lane",                required_argument,  0,                          'l'},
        {"tile",                required_argument,  0,                          't'},
        {"ncycles",             required_argument,  0,                          'n'},
        {"alternative-call",    required_argument,  0,                          'c'},
        {"index-read",          required_argument,  0,                          'i'},
        {"threep-read",         required_argument,  0,                          'e'},
        {"sample",              required_argument,  0,                          'm'},
        {"writer-command",      required_argument,  0,                          'w'},
        {"filter-control",      required_argument,  0,                          'f'},
        {"block-size",          required_argument,  0,                          'B'},
        {"demultiplex-stats",   required_argument,  0,                          'o'},
        {"color-matrix",        required_argument,  0,                          'x'},
        {"sigproc-trigger",     required_argument,  0,                          's'},
        {"minimum-polya-length",required_argument,  0,                          'M'},
        {"maximum-modifications",required_argument, 0,                          'X'},
        {"weight-T",            required_argument,  0,                          'T'},
        {"weight-ACG",          required_argument,  0,                          'A'},
        {"weight-N",            required_argument,  0,                          'N'},
        {"help",                no_argument,        0,                          'h'},
        {0, 0, 0, 0}
    };

    datadir = runid = writercmd = NULL;
    demultiplex_stats_filename = NULL;
    color_matrix_filename = NULL;
    keep_no_delimiter_flag = 0;
    lane = tile = ncycles = index_start = threep_start = -1;
    index_length = 6;
    threep_length = -1;
    barcodes = NULL;
    altcalls = NULL;
    controlinfo.name[0] = 0;
    controlinfo.first_cycle = controlinfo.read_length = -1;
    blocksize = DEFAULT_BATCH_BLOCK_SIZE;

    finderparams.max_terminal_modifications = 10;
    finderparams.min_polya_length = 5;
    finderparams.sigproc_trigger_polya_length = 10;
    memset(finderparams.weights, 0, sizeof(finderparams.weights));
    finderparams.weights[(int)'T'] = 1;
    finderparams.weights[(int)'A'] = finderparams.weights[(int)'C'] =
        finderparams.weights[(int)'G'] = -10;
    finderparams.weights[(int)'N'] = -5;

    /* XXX TEMPORARY ================ */
    rulerparams.balancer_start = 0;
    rulerparams.balancer_end = 20;
    rulerparams.balancer_minimum_occurrence = 2;
    rulerparams.balancer_num_positive_samples = 2;
    rulerparams.balancer_num_negative_samples = 4;
    rulerparams.dark_cycles_threshold = 10;
    rulerparams.max_dark_cycles = 5;
    rulerparams.polya_score_threshold = .1;
    rulerparams.downhill_extension_weight = .49;
    precalc_score_tables(&rulerparams, 20.f, .75f);
    /* XXX TEMPORARY ================ */

    while (1) {
        int option_index=0;

        c = getopt_long(argc, argv,
                        "d:r:l:t:n:s:i:e:m:w:hf:B:o:x:M:",
                        long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
            case 'd': /* --data-dir */
                datadir = strdup(optarg);
                break;

            case 'r': /* --run-id */
                runid = strdup(optarg);
                break;

            case 'l': /* --lane */
                lane = atoi(optarg);
                break;

            case 't': /* --tile */
                tile = atoi(optarg);
                break;

            case 'n': /* --ncycles */
                ncycles = atoi(optarg);
                break;

            case 'c': /* --alternative-call */
                {
#define ALTCALL_OPTION_TOKENS   2
                    struct AlternativeCallInfo *newac;
                    char *str, *saveptr;
                    char *tokens[ALTCALL_OPTION_TOKENS];
                    int j;

                    saveptr = NULL;
                    for (j = 0, str = optarg; j < ALTCALL_OPTION_TOKENS; j++, str = NULL) {
                        tokens[j] = strtok_r(str, ",", &saveptr);
                        if (tokens[j] == NULL)
                            break;
                    }

                    if (j != ALTCALL_OPTION_TOKENS) {
                        fprintf(stderr, "An alternative call specified with illegal format.\n");
                        return -1;
                    }

                    newac = malloc(sizeof(struct AlternativeCallInfo));
                    if (newac == NULL) {
                        perror("main");
                        return -1;
                    }

                    newac->filename = strdup(tokens[0]);
                    newac->first_cycle = atoi(tokens[1]) - 1;
                    newac->reader = NULL;

                    if (newac->filename == NULL) {
                        perror("main");
                        return -1;
                    }

                    newac->next = altcalls;
                    altcalls = newac;
                }
                break;

            case 'i': /* --index-read */
                {
#define INDEX_OPTION_TOKENS   2
                    char *str, *saveptr;
                    char *tokens[INDEX_OPTION_TOKENS];
                    int j;

                    saveptr = NULL;
                    for (j = 0, str = optarg; j < INDEX_OPTION_TOKENS; j++, str = NULL) {
                        tokens[j] = strtok_r(str, ",", &saveptr);
                        if (tokens[j] == NULL)
                            break;
                    }

                    if (j != INDEX_OPTION_TOKENS) {
                        fprintf(stderr, "Index read range specified in illegal format.\n");
                        return -1;
                    }

                    index_start = atoi(tokens[0]) - 1;
                    index_length = atoi(tokens[1]);
                }
                break;

            case 'e': /* --threep-read */
                {
#define THREEP_OPTION_TOKENS   2
                    char *str, *saveptr;
                    char *tokens[THREEP_OPTION_TOKENS];
                    int j;

                    saveptr = NULL;
                    for (j = 0, str = optarg; j < THREEP_OPTION_TOKENS; j++, str = NULL) {
                        tokens[j] = strtok_r(str, ",", &saveptr);
                        if (tokens[j] == NULL)
                            break;
                    }

                    if (j != THREEP_OPTION_TOKENS) {
                        fprintf(stderr, "3' read range specified in illegal format.\n");
                        return -1;
                    }

                    threep_start = atoi(tokens[0]) - 1;
                    threep_length = atoi(tokens[1]);
                }
                break;

            case 'm': /* --sample */
                {
#define SAMPLE_OPTION_TOKENS        8
                    struct BarcodeInfo *newbc;
                    char *str;
                    char *tokens[SAMPLE_OPTION_TOKENS];
                    int j;

                    for (j = 0, str = optarg; j < SAMPLE_OPTION_TOKENS; j++) {
                        tokens[j] = strsep(&str, ",");
                        if (tokens[j] == NULL)
                            break;
                    }

                    if (j != SAMPLE_OPTION_TOKENS) {
                        fprintf(stderr, "j=%d\n", j);
                        fprintf(stderr, "A sample specified with illegal format.\n");
                        return -1;
                    }

                    newbc = malloc(sizeof(struct BarcodeInfo));
                    if (newbc == NULL) {
                        perror("main");
                        return -1;
                    }

                    newbc->name = strdup(tokens[0]);
                    newbc->index = strdup(tokens[1]);
                    newbc->maximum_index_mismatches = atoi(tokens[2]);
                    newbc->delimiter = strdup(tokens[3]);
                    newbc->delimiter_pos = atoi(tokens[4]) - 1;
                    newbc->delimiter_length = strlen(newbc->delimiter);
                    newbc->maximum_delimiter_mismatches = atoi(tokens[5]);
                    newbc->fingerprint = strdup(tokens[6]);
                    newbc->fingerprint_length = strlen(newbc->fingerprint);
                    newbc->maximum_fingerprint_mismatches = atoi(tokens[7]);
                    newbc->clusters_mm0 = newbc->clusters_mm1 = 0;
                    newbc->clusters_mm2plus = newbc->clusters_nodelim = 0;
                    newbc->clusters_fpmismatch = 0;
                    newbc->stream = NULL;

                    if (newbc->name == NULL || newbc->index == NULL ||
                                newbc->delimiter == NULL) {
                        perror("main");
                        return -1;
                    }

                    if (newbc->maximum_index_mismatches < 0 ||
                            newbc->maximum_index_mismatches >= strlen(newbc->index)) {
                        fprintf(stderr, "Maximum index mismatches were out of range for "
                                        "sample `%s'.\n", newbc->name);
                        return -1;
                    }

                    if (newbc->delimiter_length > 0 && (
                            newbc->maximum_delimiter_mismatches < 0 ||
                            newbc->maximum_delimiter_mismatches >= strlen(newbc->delimiter))) {
                        fprintf(stderr, "Maximum delimiter mismatches were out of range "
                                        "for sample `%s'.\n", newbc->name);
                        return -1;
                    }

                    if (newbc->fingerprint_length > 0 && (
                            newbc->maximum_fingerprint_mismatches < 0 ||
                            newbc->maximum_fingerprint_mismatches >= newbc->fingerprint_length)) {
                        fprintf(stderr, "Maximum fingerprint mismatches were out of range "
                                        "for sample `%s'.\n", newbc->name);
                        return -1;
                    }

                    newbc->next = barcodes;
                    barcodes = newbc;
                }
                break;

            case 'w': /* --writer-command */
                {
                    const char *replptr;

                    if (strchr(optarg, '%') != NULL) {
                        fprintf(stderr, "The writer command must not include the `%%' "
                                        "character.\n");
                        return -1;
                    }

                    replptr = strstr(optarg, "XX");
                    if (replptr == NULL) {
                        fprintf(stderr, "The writer command must include an incidence of `XX' "
                                        "to be replaced with sample names.\n");
                        return -1;
                    }

                    writercmd = strdup(optarg);
                    writercmd[replptr - optarg] = '%';
                    writercmd[replptr - optarg + 1] = 's';
                }
                break;

            case 'f':
                {
#define FILTER_CONTROL_OPTION_TOKENS    3
                    char *str, *saveptr;
                    char *tokens[FILTER_CONTROL_OPTION_TOKENS];
                    int j;

                    saveptr = NULL;
                    for (j = 0, str = optarg; j < FILTER_CONTROL_OPTION_TOKENS;
                                              j++, str = NULL) {
                        tokens[j] = strtok_r(str, ",", &saveptr);
                        if (tokens[j] == NULL)
                            break;
                    }

                    if (j != FILTER_CONTROL_OPTION_TOKENS) {
                        fprintf(stderr, "A sample specified with illegal format.\n");
                        return -1;
                    }

                    strncpy(controlinfo.name, tokens[0], CONTROL_NAME_MAX);
                    controlinfo.first_cycle = atoi(tokens[1]) - 1;
                    controlinfo.read_length = atoi(tokens[2]);
                    controlinfo.barcode = NULL;
                }
                break;

            case 'B': /* --block-size */
                blocksize = atoi(optarg);
                break;

            case 'o': /* --demultiplex-stats */
                demultiplex_stats_filename = strdup(optarg);
                break;

            case 'x': /* --color-matrix */
                color_matrix_filename = strdup(optarg);
                break;

            case 's': /* --sigproc-trigger */
                finderparams.sigproc_trigger_polya_length = atoi(optarg);
                break;

            case 'M': /* --minimum-polya-length */
                finderparams.min_polya_length = atoi(optarg);
                break;

            case 'X': /* --maximum-modifications */
                finderparams.max_terminal_modifications = atoi(optarg);
                break;

            case 'T': /* --weight-T */
                finderparams.weights[(int)'T'] = atoi(optarg);
                break;

            case 'A': /* --weight-ACG */
                finderparams.weights[(int)'A'] =
                    finderparams.weights[(int)'C'] =
                    finderparams.weights[(int)'G'] = atoi(optarg);
                break;

            case 'N': /* --weight-N */
                finderparams.weights[(int)'N'] = atoi(optarg);
                break;

            case 'h':
                usage(argv[0]);
                exit(0);
                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort();
        }
    }

    if (datadir == NULL) {
        usage(argv[0]);
        fprintf(stderr, "--data-dir is not set.\n");
        return -1;
    }

    if (runid == NULL) {
        usage(argv[0]);
        fprintf(stderr, "--run-id is not set.\n");
        return -1;
    }

    if (lane < 0) {
        usage(argv[0]);
        fprintf(stderr, "--lane is not set.\n");
        return -1;
    }

    if (tile < 0) {
        usage(argv[0]);
        fprintf(stderr, "--tile is not set.\n");
        return -1;
    }

    if (ncycles < 0) {
        usage(argv[0]);
        fprintf(stderr, "--ncycles is not set.\n");
        return -1;
    }

    if (index_start < 0) {
        usage(argv[0]);
        fprintf(stderr, "--index-read is not set.\n");
        return -1;
    }

    if (threep_start < 0) {
        usage(argv[0]);
        fprintf(stderr, "--threep-read is not set.\n");
        return -1;
    }

    if (barcodes == NULL) {
        usage(argv[0]);
        fprintf(stderr, "--sample is not set.\n");
        return -1;
    }

    if (writercmd == NULL) {
        usage(argv[0]);
        fprintf(stderr, "--writer-command is not set.\n");
        return -1;
    }

    if (color_matrix_filename == NULL) {
        usage(argv[0]);
        fprintf(stderr, "--color-matrix is not set.\n");
        return -1;
    }

    if (load_color_matrix(rulerparams.colormatrix,
                          color_matrix_filename) < 0)
        return -1;

    if (controlinfo.first_cycle >= 0) {
        struct BarcodeInfo *control;

        control = malloc(sizeof(struct BarcodeInfo));
        memset(control, 0, sizeof(struct BarcodeInfo));
        control->name = strdup(controlinfo.name);

        control->index = malloc(index_length + 1);
        memset(control->index, 'X', index_length);

        control->delimiter = strdup("");
        control->delimiter_pos = -1;
        control->maximum_index_mismatches = index_length;
        control->maximum_delimiter_mismatches = -1;
        control->fingerprint = strdup("");
        control->next = barcodes;
        controlinfo.barcode = control;

        barcodes = control;
    }

    /* Add a fallback barcode entry for "Unknown" samples. */
    {
        struct BarcodeInfo *fallback;

        fallback = malloc(sizeof(struct BarcodeInfo));
        memset(fallback, 0, sizeof(struct BarcodeInfo));
        fallback->name = strdup("Unknown");

        fallback->index = malloc(index_length + 1);
        memset(fallback->index, 'X', index_length);

        fallback->delimiter = strdup("");
        fallback->delimiter_pos = -1;
        fallback->maximum_index_mismatches = index_length;
        fallback->maximum_delimiter_mismatches = -1;
        fallback->fingerprint = strdup("");
        fallback->next = barcodes;

        barcodes = fallback;
    }

    {
        struct ControlFilterInfo *pctlinfo;

        if (controlinfo.first_cycle >= 0)
            pctlinfo = &controlinfo;
        else
            pctlinfo = NULL;

        r = process(datadir, runid, lane, tile, ncycles,
                    threep_start, threep_length, index_start, index_length,
                    barcodes, writercmd, altcalls,
                    pctlinfo, &finderparams, &rulerparams, blocksize,
                    keep_no_delimiter_flag);
    }
    
    free(datadir);
    free(runid);
    free(writercmd);

    if (r == 0 && demultiplex_stats_filename != NULL) {
        r = write_demultiplexing_statistics(demultiplex_stats_filename, barcodes);
        free(demultiplex_stats_filename);
    }

    while (barcodes != NULL) {
        struct BarcodeInfo *bk;

        if (barcodes->stream != NULL)
            pclose(barcodes->stream);

        free(barcodes->name);
        free(barcodes->index);
        free(barcodes->delimiter);
        bk = barcodes->next;
        free(barcodes);
        barcodes = bk;
    }

    while (altcalls != NULL) {
        struct AlternativeCallInfo *ac;

        ac = altcalls->next;
        free(altcalls);
        altcalls = ac;
    }

    return r;
}
