/*
 * tailseq-retrieve-signals.c
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
#include "tailseq-retrieve-signals.h"


#define DEFAULT_BATCH_BLOCK_SIZE    262144
/* BCL - 2 bytes per cluster per cycle
 * CIF - 8 bytes per cluster per cycle
 * 3080 bytes per cluster for regular 51+251 set-up.
 * then, 262144 clusters takes 770 MiB up for storing data. */


static struct BarcodeInfo *
assign_barcode(const char *indexseq, int barcode_length, struct BarcodeInfo *barcodes)
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

    if (bestidx == NULL || secondbestfound || bestmismatches > bestidx->maximum_mismatches)
        return NULL;
    else
        return bestidx;
}


static int
find_delimiter_end_position(const char *sequence, struct BarcodeInfo *barcode)
{
    const char *delimpos;

    delimpos = sequence + barcode->delimiter_pos;

    if (memcmp(delimpos, barcode->delimiter, barcode->delimiter_length) == 0)
        return barcode->delimiter_pos + barcode->delimiter_length;

    if (memcmp(delimpos - 1, barcode->delimiter, barcode->delimiter_length) == 0)
        return barcode->delimiter_pos + barcode->delimiter_length - 1;

    if (memcmp(delimpos + 1, barcode->delimiter, barcode->delimiter_length) == 0)
        return barcode->delimiter_pos + barcode->delimiter_length + 1;

    return -1;
}


static int
demultiplex_and_write(const char *laneid, int tile, int ncycles, uint32_t firstclusterno,
                      int scalefactor, int barcode_start, int barcode_length,
                      struct BarcodeInfo *barcodes,
                      struct CIFData **intensities, struct BCLData **basecalls,
                      struct ControlFilterInfo *control_info, int keep_no_delimiter)
{
    uint32_t cycleno, clusterno;
    uint32_t clustersinblock;
    char sequence_formatted[ncycles+1], quality_formatted[ncycles+1];
    char intensity_formatted[ncycles*8+1];
    struct BarcodeInfo *noncontrol_barcodes;

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

    for (clusterno = 0; clusterno < clustersinblock; clusterno++) {
        struct BarcodeInfo *bc;
        int delimiter_end;

        format_basecalls(sequence_formatted, quality_formatted, basecalls, ncycles, clusterno);
        format_intensity(intensity_formatted, intensities, ncycles, clusterno, scalefactor);

        bc = assign_barcode(sequence_formatted + barcode_start, barcode_length,
                            noncontrol_barcodes);
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

        delimiter_end = find_delimiter_end_position(sequence_formatted, bc);
        if (bc->delimiter_length > 0 && delimiter_end < 0 && !keep_no_delimiter)
            continue;

        if (fprintf(bc->stream, "%s%04d\t%d\t%d\t%s\t%s\t%s\n", laneid, tile,
                    firstclusterno + clusterno, delimiter_end, sequence_formatted,
                    quality_formatted, intensity_formatted) < 0) {
            perror("demultiplex_and_write");

            if (control_seq != NULL)
                free(control_seq);
            return -1;
        }
    }

    if (control_seq != NULL)
        free(control_seq);

    return 0;
}


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
        int scalefactor, int barcode_start, int barcode_length,
        struct BarcodeInfo *barcodes, const char *writercmd,
        struct AlternativeCallInfo *altcalls, struct ControlFilterInfo *control_info,
        int blocksize, int keep_no_delimiter)
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
        if (demultiplex_and_write(laneid, tile, ncycles, nclusters - clusters_to_go,
                                  scalefactor, barcode_start, barcode_length, barcodes,
                                  intensities, basecalls, control_info,
                                  keep_no_delimiter) == -1)
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


static void
usage(const char *prog)
{
    printf("\
tailseq-retrieve-signals 2.0\
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
\n  -b,  --barcode-start=NUM          the first cycle of index read.\
\n  -a,  --barcode-length=NUM         length of index read.\
\n  -w,  --writer-command=COMMAND     shell command to run to write .sqi\
\n                                    output (usually a stream compressor.)\
\n\
\nOptional parameters:\
\n  -s,  --signal-scale=NUM           number of digits in radix 2 to scale\
\n                                    signal intensity down. (default: 0)\
\n  -c,  --alternative-call=SPEC      FASTQ file to replace base calls.\
\n                                    specify as \"filename,first_cycle\".\
\n  -m,  --sample=SPEC                sample information as formatted in\
\n            \"name,index,maximum_mismatches,delimiter,delimiter_position\".\
\n       --keep-no-delimiter          don't skip reads where a delimiter\
\n                                    is not found.\
\n  -f,  --filter-control=SPEC        sort out PhiX control reads by sequence\
\n                                    in \"name,first_cycle,length\".\
\n  -B,  --block-size=NUM             number of clusters in a processing unit\
\n                                    (default: 262144).\
\n\
\nAll cycle numbers or positions are in 1-based inclusive system.\
\n\
\nMail bug reports and suggestions to Hyeshik Chang <hyeshik@snu.ac.kr>.\n", prog);
}


int
main(int argc, char *argv[])
{
    char *datadir, *runid, *writercmd;
    int keep_no_delimiter_flag;
    int lane, tile, ncycles, blocksize;
    int barcode_start, barcode_length;
    int scalefactor;
    struct BarcodeInfo *barcodes;
    struct AlternativeCallInfo *altcalls;
    struct ControlFilterInfo controlinfo;
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
        {"signal-scale",        required_argument,  0,                          's'},
        {"barcode-start",       required_argument,  0,                          'b'},
        {"barcode-length",      required_argument,  0,                          'a'},
        {"sample",              required_argument,  0,                          'm'},
        {"writer-command",      required_argument,  0,                          'w'},
        {"filter-control",      required_argument,  0,                          'f'},
        {"block-size",          required_argument,  0,                          'B'},
        {"help",                no_argument,        0,                          'h'},
        {0, 0, 0, 0}
    };

    datadir = runid = writercmd = NULL;
    keep_no_delimiter_flag = 0;
    lane = tile = ncycles = barcode_start = -1;
    scalefactor = 0;
    barcode_length = 6;
    barcodes = NULL;
    altcalls = NULL;
    controlinfo.name[0] = 0;
    controlinfo.first_cycle = controlinfo.read_length = -1;
    blocksize = DEFAULT_BATCH_BLOCK_SIZE;

    while (1) {
        int option_index=0;

        c = getopt_long(argc, argv, "d:r:l:t:n:s:b:a:m:w:hf:B:", long_options, &option_index);

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

            case 's': /* --signal-scale */
                scalefactor = atoi(optarg);
                break;

            case 'b': /* --barcode-start */
                barcode_start = atoi(optarg) - 1;
                break;

            case 'a': /* --barcode-length */
                barcode_length = atoi(optarg);
                break;

            case 'm': /* --sample */
                {
#define SAMPLE_OPTION_TOKENS        5
                    struct BarcodeInfo *newbc;
                    char *str, *saveptr;
                    char *tokens[SAMPLE_OPTION_TOKENS];
                    int j;

                    for (j = 0, str = optarg; j < SAMPLE_OPTION_TOKENS; j++, str = NULL) {
                        tokens[j] = strtok_r(str, ",", &saveptr);
                        if (tokens[j] == NULL)
                            break;
                    }

                    if (j != SAMPLE_OPTION_TOKENS) {
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
                    newbc->maximum_mismatches = atoi(tokens[2]);
                    newbc->delimiter = strdup(tokens[3]);
                    newbc->delimiter_pos = atoi(tokens[4]) - 1;
                    newbc->delimiter_length = strlen(newbc->delimiter);
                    newbc->stream = NULL;

                    if (newbc->name == NULL || newbc->index == NULL ||
                                newbc->delimiter == NULL) {
                        perror("main");
                        return -1;
                    }

                    if (newbc->maximum_mismatches < 0 ||
                            newbc->maximum_mismatches >= strlen(newbc->index)) {
                        fprintf(stderr, "Maximum mismatches were out of range for sample "
                                        "`%s'.\n", newbc->name);
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

    if (barcode_start < 0) {
        usage(argv[0]);
        fprintf(stderr, "--barcode-start is not set.\n");
        return -1;
    }

    if (barcode_length < 0) {
        usage(argv[0]);
        fprintf(stderr, "--barcode-length is not set.\n");
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

    if (controlinfo.first_cycle >= 0) {
        struct BarcodeInfo *control;

        control = malloc(sizeof(struct BarcodeInfo));
        control->name = strdup(controlinfo.name);

        control->index = malloc(barcode_length + 1);
        memset(control->index, 'X', barcode_length);
        control->index[barcode_length] = '\0';

        control->delimiter = strdup("");
        control->delimiter_pos = -1;
        control->delimiter_length = 0;
        control->maximum_mismatches = barcode_length;
        control->stream = NULL;
        control->next = barcodes;
        controlinfo.barcode = control;

        barcodes = control;
    }

    /* Add a fallback barcode entry for "Unknown" samples. */
    {
        struct BarcodeInfo *fallback;

        fallback = malloc(sizeof(struct BarcodeInfo));
        fallback->name = strdup("Unknown");

        fallback->index = malloc(barcode_length + 1);
        memset(fallback->index, 'X', barcode_length);
        fallback->index[barcode_length] = '\0';

        fallback->delimiter = strdup("");
        fallback->delimiter_pos = -1;
        fallback->delimiter_length = 0;
        fallback->maximum_mismatches = barcode_length;
        fallback->stream = NULL;
        fallback->next = barcodes;

        barcodes = fallback;
    }

    {
        struct ControlFilterInfo *pctlinfo;

        if (controlinfo.first_cycle >= 0)
            pctlinfo = &controlinfo;
        else
            pctlinfo = NULL;

        r = process(datadir, runid, lane, tile, ncycles, scalefactor, barcode_start,
                    barcode_length, barcodes, writercmd, altcalls,
                    pctlinfo, blocksize, keep_no_delimiter_flag);
    }
    
    free(datadir);
    free(runid);
    free(writercmd);

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
