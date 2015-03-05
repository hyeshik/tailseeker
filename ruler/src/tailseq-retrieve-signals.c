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
#include <endian.h>
#include <zlib.h>

#define NUM_CHANNELS        4
#define NOCALL_QUALITY      2
#define NOCALL_BASE         'N'
#define PHRED_BASE          33
#define MAX_INDEX_LENGTH    16


struct IntensitySet {
    uint16_t value[NUM_CHANNELS];
};

struct CIFData {
    size_t nclusters;
    struct IntensitySet intensity[1];
};

struct CIFHeader {
    char filemagic[3];
    uint8_t version;
    uint8_t datasize;
    uint16_t first_cycle;
    uint16_t ncycles;
    uint32_t nclusters;
};

struct BCLData {
    size_t nclusters;
    uint8_t *base;
    uint8_t *quality;
};

struct BarcodeInfo;
struct BarcodeInfo {
    char *name;
    char *index;
    char *delimiter;
    int delimiter_pos;
    int delimiter_length;
    int maximum_mismatches;
    FILE *stream;
    struct BarcodeInfo *next;
};

struct AlternativeCallInfo;
struct AlternativeCallInfo {
    char *filename;
    int first_cycle; /* the last cycle number will be determined from the file content */
    struct AlternativeCallInfo *next;
};

static const char CALL_BASES[4] = "ACGT";
static const char BASE64_ENCODE_TABLE[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef"
                                          "ghijklmnopqrstuvwxyz0123456789+/";


static struct CIFData *
load_cif_file(const char *filename)
{
    FILE *fp;
    struct CIFHeader header;
    struct CIFData *data;

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "load_cif_file: Can't open file %s.\n", filename);
        return NULL;
    }

    if (fread(header.filemagic, 3, 1, fp) < 1)
        goto onError;

    if (memcmp(header.filemagic, "CIF", 3) != 0) {
        fprintf(stderr, "File %s does not seem to be a CIF file.", filename);
        goto onError;
    }

    /* read version and data size */
    if (fread(&header.version, 2, 1, fp) < 1) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        goto onError;
    }
    if (header.version != 1) {
        fprintf(stderr, "Unsupported CIF version: %d.\n", header.version);
        goto onError;
    }
    if (header.datasize != 2) {
        fprintf(stderr, "Unsupported data size (%d).\n", header.datasize);
        goto onError;
    }

    /* read first cycle and number of cycles */
    if (fread(&header.first_cycle, 2, 2, fp) < 2) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        goto onError;
    }
    if (fread(&header.nclusters, 4, 1, fp) < 1) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        goto onError;
    }

    header.first_cycle = le16toh(header.first_cycle);
    header.ncycles = le16toh(header.ncycles);
    header.nclusters = le32toh(header.nclusters);

    data = malloc(sizeof(int) + sizeof(struct IntensitySet) * header.nclusters);
    if (data == NULL) {
        perror("load_cif_file");
        goto onError;
    }

    data->nclusters = header.nclusters;
    if (fread(data->intensity, sizeof(struct IntensitySet), header.nclusters, fp) !=
            header.nclusters) {
        fprintf(stderr, "Not all data were loaded from %s.", filename);
        goto onError;
    }

    fclose(fp);

    return data;

  onError:
    fclose(fp);
    return NULL;
}


static void
free_cif_data(struct CIFData *data)
{
    free(data);
}


static struct BCLData *
load_bcl_file(const char *filename)
{
    struct BCLData *bcl=NULL;
    uint32_t nclusters, i;
    uint8_t *base, *quality;
    FILE *fp;

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "load_bcl_file: Can't open file %s.\n", filename);
        return NULL;
    }

    if (fread(&nclusters, sizeof(nclusters), 1, fp) < 1) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        goto onError;
    }

    nclusters = le32toh(nclusters);
    bcl = malloc(sizeof(struct BCLData));
    if (bcl == NULL) {
        perror("load_bcl_file");
        goto onError;
    }

    bcl->nclusters = nclusters;
    base = bcl->base = malloc(nclusters);
    quality = bcl->quality = malloc(nclusters);
    if (base == NULL || quality == NULL) {
        perror("load_bcl_file");
        goto onError;
    }

    if (fread(base, nclusters, 1, fp) < 1) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        printf("nclusters=%u file=%s\n", nclusters, filename);
        goto onError;
    }

    fclose(fp);

    for (i = 0; i < nclusters; i++, base++, quality++) {
        if (*base == 0) {
            *quality = NOCALL_QUALITY + PHRED_BASE;
            *base = NOCALL_BASE;
        }
        else {
            *quality = (*base >> 2) + PHRED_BASE;
            *base = CALL_BASES[*base & 3];
        }
    }

    return bcl;

  onError:
    if (bcl != NULL) {
        if (bcl->base != NULL)
            free(bcl->base);
        if (bcl->quality != NULL)
            free(bcl->quality);
        free(bcl);
    }

    fclose(fp);

    return NULL;
}


static void
free_bcl_data(struct BCLData *data)
{
    free(data->base);
    free(data->quality);
    free(data);
}


static int
load_alternative_basecalls(const char *filename, struct BCLData **basecalls,
                           uint32_t nclusters)
{
    char buf[BUFSIZ];
    gzFile hdl;
    uint32_t clusterno;
    int ncycles, j;

    ncycles = -1;

    hdl = gzopen(filename, "rt");
    if (hdl == NULL) {
        perror("load_alternative_basecalls");
        return -1;
    }

    for (clusterno = 0; clusterno < nclusters; clusterno++) {
        size_t seqlen;

        /* header 1 */
        if (gzgets(hdl, buf, BUFSIZ) == NULL) {
            fprintf(stderr, "Not enough sequences from %s\n", filename);
            return -1;
        }

        if (buf[0] != '@') {
            fprintf(stderr, "%s is not in FASTQ format.\n", filename);
            return -1;
        }

        /* sequence */
        if (gzgets(hdl, buf, BUFSIZ) == NULL) {
            fprintf(stderr, "Unexpected EOF while reading %s\n", filename);
            return -1;
        }

        seqlen = strlen(buf) - 1;
        if (ncycles == -1) {
            /* Initialize here */
            ncycles = seqlen;

            for (j = 0; j < ncycles; j++) {
                basecalls[j] = malloc(sizeof(struct BCLData));
                if (basecalls[j] == NULL) {
                    perror("load_alternative_basecalls");
                    gzclose(hdl);
                    basecalls[j]->base = basecalls[j]->quality = NULL;
                    return -1;
                }

                basecalls[j]->nclusters = nclusters;
                basecalls[j]->base = malloc(nclusters);
                basecalls[j]->quality = malloc(nclusters);
                if (basecalls[j]->base == NULL || basecalls[j]->quality == NULL) {
                    perror("load_alternative_basecalls");
                    gzclose(hdl);
                    return -1;
                }
            }
        }
        else if (seqlen != ncycles) {
            fprintf(stderr, "Sequence length in the FASTQ is inconsistent.\n");
            return -1;
        }

        for (j = 0; j < ncycles; j++)
            basecalls[j]->base[clusterno] = buf[j];

        /* header 2 */
        if (gzgets(hdl, buf, BUFSIZ) == NULL) {
            fprintf(stderr, "Unexpected EOF while reading %s\n", filename);
            return -1;
        }

        if (buf[0] != '+') {
            fprintf(stderr, "%s is not in FASTQ format.\n", filename);
            return -1;
        }

        /* quality */
        if (gzgets(hdl, buf, BUFSIZ) == NULL) {
            fprintf(stderr, "Unexpected EOF while reading %s\n", filename);
            return -1;
        }

        seqlen = strlen(buf) - 1;
        if (seqlen != ncycles) {
            fprintf(stderr, "Sequence length in the FASTQ is inconsistent.\n");
            return -1;
        }

        for (j = 0; j < ncycles; j++)
            basecalls[j]->quality[clusterno] = buf[j];
    }

    if (gzgets(hdl, buf, BUFSIZ) != NULL) {
        fprintf(stderr, "Extra sequences found in %s\n", filename);
        return -1;
    }

    gzclose(hdl);

    return ncycles;
}


static int
load_intensities_and_basecalls(const char *datadir, const char *laneid, int lane, int tile,
                               int ncycles, struct AlternativeCallInfo *altcalls,
                               struct CIFData **intensities, struct BCLData **basecalls)
{
    char path[PATH_MAX];
    char altcall_loaded[ncycles];
    uint32_t cycleno;
    struct AlternativeCallInfo *paltcall;
    int bcllaststart;

    memset(altcall_loaded, 0, ncycles);
    bcllaststart = -1;

    printf("[%s%d] - CIF data source: %s\n", laneid, tile, datadir);

    for (cycleno = 0; cycleno < ncycles; cycleno++) {
        snprintf(path, PATH_MAX, "%s/L%03d/C%d.1/s_%d_%04d.cif", datadir, lane, cycleno+1,
                 lane, tile);

        intensities[cycleno] = load_cif_file(path);
        if (intensities[cycleno] == NULL)
            return -1;

        if (altcall_loaded[cycleno])
            continue;

        for (paltcall = altcalls; paltcall != NULL; paltcall = paltcall->next) {
            if (paltcall->first_cycle == cycleno) {
                int loaded, j;

                loaded = load_alternative_basecalls(paltcall->filename, basecalls + cycleno,
                                                    intensities[cycleno]->nclusters);
                if (loaded < 0)
                    return -1;

                for (j = cycleno; j < cycleno + loaded; j++)
                    altcall_loaded[j] = 1;

                if (bcllaststart >= 0) {
                    printf("[%s%d] - BCL basecall source [%d-%d]: %s\n", laneid, tile,
                            bcllaststart+1, cycleno, datadir);
                    bcllaststart = -1;
                }

                printf("[%s%d] - FASTQ basecall source [%d-%d]: %s\n", laneid, tile,
                        cycleno+1, cycleno+loaded, paltcall->filename);

                break;
            }
        }

        if (paltcall != NULL) /* exited after loading an alternative basecall file. */
            continue;

        snprintf(path, PATH_MAX, "%s/BaseCalls/L%03d/C%d.1/s_%d_%04d.bcl", datadir, lane,
                 cycleno+1, lane, tile);

        basecalls[cycleno] = load_bcl_file(path);
        if (basecalls[cycleno] == NULL)
            return -1;

        if (bcllaststart < 0)
            bcllaststart = cycleno;
    }

    if (bcllaststart >= 0)
        printf("[%s%d] - BCL basecall source [%d-%d]: %s\n", laneid, tile,
                bcllaststart+1, cycleno, datadir);

    return 0;
}


static void
format_basecalls(char *seq, char *qual, struct BCLData **basecalls,
                 int ncycles, uint32_t clusterno)
{
    uint32_t i;

    for (i = 0; i < ncycles; i++) {
        *seq++ = basecalls[i]->base[clusterno];
        *qual++ = basecalls[i]->quality[clusterno];
    }

    *seq = *qual = 0;
}


static void
format_intensity(char *inten, struct CIFData **intensities,
                 int ncycles, uint32_t clusterno, double scalefactor)
{
    uint32_t i;
    int value, chan;

    for (i = 0; i < ncycles; i++) {
        for (chan = 0; chan < NUM_CHANNELS; chan++) {
            value = scalefactor * intensities[i]->intensity[clusterno].value[chan];
            value += 255;
            if (value >= 4096)
                value = 4095;
            else if (value < 0)
                value = 0;

            *inten++ = BASE64_ENCODE_TABLE[value >> 6];
            *inten++ = BASE64_ENCODE_TABLE[value & 63];
        }
    }

    *inten = 0;
}


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
        return barcodes; /* which is "Unknown". */
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
demultiplex_and_write(const char *laneid, int tile, int ncycles, double scalefactor,
                      int barcode_start, int barcode_length,
                      struct BarcodeInfo *barcodes,
                      struct CIFData **intensities, struct BCLData **basecalls,
                      int keep_no_delimiter)
{
    uint32_t cycleno, clusterno;
    uint32_t totalclusters;
    char sequence_formatted[ncycles+1], quality_formatted[ncycles+1];
    char intensity_formatted[ncycles*8+1];

    totalclusters = intensities[0]->nclusters;
    for (cycleno = 0; cycleno < ncycles; cycleno++)
        if (totalclusters != intensities[cycleno]->nclusters ||
                totalclusters != basecalls[cycleno]->nclusters) {
            fprintf(stderr, "Inconsistent number of clusters in cycle %d.\n", cycleno);
            return -1;
        }

    for (clusterno = 0; clusterno < totalclusters; clusterno++) {
        struct BarcodeInfo *bc;
        int delimiter_end;

        format_basecalls(sequence_formatted, quality_formatted, basecalls, ncycles, clusterno);
        format_intensity(intensity_formatted, intensities, ncycles, clusterno, scalefactor);

        bc = assign_barcode(sequence_formatted + barcode_start, barcode_length, barcodes);
        delimiter_end = find_delimiter_end_position(sequence_formatted, bc);
        if (bc->delimiter_length > 0 && delimiter_end < 0 && !keep_no_delimiter)
            continue;

        if (fprintf(bc->stream, "%s%04d\t%d\t0\t%s\t%s\t%s\t%d\n", laneid, tile, clusterno,
             sequence_formatted, quality_formatted, intensity_formatted, delimiter_end) < 0) {
            perror("demultiplex_and_write");
            return -1;
        }
    }

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
        if (barcodes->stream != NULL)
            pclose(barcodes->stream);
}


static int
process(const char *datadir, const char *laneid, int lane, int tile, int ncycles,
        double scalefactor, int barcode_start, int barcode_length,
        struct BarcodeInfo *barcodes, const char *writercmd,
        struct AlternativeCallInfo *altcalls, int keep_no_delimiter)
{
    struct CIFData **intensities;
    struct BCLData **basecalls;
    int cycleno;

    intensities = malloc(sizeof(struct CIFData *) * ncycles);
    if (intensities == NULL) {
        perror("process");
        return -1;
    }

    basecalls = malloc(sizeof(struct BCLData *) * ncycles);
    if (basecalls == NULL) {
        free(intensities);
        perror("process");
        return -1;
    }

    memset(intensities, 0, sizeof(struct CIFData *) * ncycles);
    memset(basecalls, 0, sizeof(struct BCLData *) * ncycles);

    printf("[%s%d] Loading CIF and BCL files\n", laneid, tile);

    if (load_intensities_and_basecalls(datadir, laneid, lane, tile, ncycles, altcalls,
                                       intensities, basecalls) == -1)
        goto onError;

    printf("[%s%d] Opening writer subprocesses\n", laneid, tile);
    if (spawn_writers(writercmd, barcodes) == -1)
        goto onError;

    printf("[%s%d] Demultiplexing and writing\n", laneid, tile);

    if (demultiplex_and_write(laneid, tile, ncycles, scalefactor, barcode_start,
                              barcode_length, barcodes, intensities, basecalls,
                              keep_no_delimiter) == -1)
        goto onError;

    for (cycleno = 0; cycleno < ncycles; cycleno++) {
        free_cif_data(intensities[cycleno]);
        free_bcl_data(basecalls[cycleno]);
    }

    terminate_writers(barcodes);

    printf("[%s%d] Finished.\n", laneid, tile);

    return 0;

  onError:
    for (cycleno = 0; cycleno < ncycles; cycleno++) {
        if (intensities[cycleno] != NULL)
            free_cif_data(intensities[cycleno]);

        if (basecalls[cycleno] != NULL)
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
tailseq-retrieve-signals 1.0\
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
\n  -s,  --signal-scale=NUM           coefficient to multiply to signal\
\n                                    intensity values.\
\n  -b,  --barcode-start=NUM          the first cycle of index read.\
\n  -a,  --barcode-length=NUM         length of index read.\
\n  -w,  --writer-command=COMMAND     shell command to run to write .sqi\
\n                                    output (usually a stream compressor.)\
\n\
\nOptional parameters:\
\n  -c,  --alternative-call=SPEC      FASTQ file to replace base calls.\
\n                                    specify as \"filename,first_cycle\".\
\n  -m,  --sample=SPEC                sample information as formatted in\
\n            \"name,index,maximum_mismatches,delimiter,delimiter_position\".\
\n       --keep-no-delimiter          don't skip reads where a delimiter\
\n                                    is not found.\
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
    int lane, tile, ncycles;
    int barcode_start, barcode_length;
    double scalefactor;
    struct BarcodeInfo *barcodes;
    struct AlternativeCallInfo *altcalls;
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
        {"help",                no_argument,        0,                          'h'},
        {0, 0, 0, 0}
    };

    datadir = runid = writercmd = NULL;
    keep_no_delimiter_flag = 0;
    lane = tile = ncycles = barcode_start = -1;
    scalefactor = -1.0;
    barcode_length = 6;
    barcodes = NULL;
    altcalls = NULL;

    while (1) {
        int option_index=0;

        c = getopt_long(argc, argv, "d:r:l:t:n:s:b:a:m:w:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
            case 0:
                printf("zero??\n");
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf("option %s", long_options[option_index].name);
                if (optarg)
                    printf(" with arg %s", optarg);
                printf("\n");
                break;

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

                    if (newac->filename == NULL) {
                        perror("main");
                        return -1;
                    }

                    newac->next = altcalls;
                    altcalls = newac;
                }
                break;

            case 's': /* --signal-scale */
                scalefactor = atof(optarg);
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

    if (scalefactor < 0.) {
        usage(argv[0]);
        fprintf(stderr, "--signal-scale is not set.\n");
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

    r = process(datadir, runid, lane, tile, ncycles, scalefactor, barcode_start,
                barcode_length, barcodes, writercmd, altcalls, keep_no_delimiter_flag);
    
    free(datadir);
    free(runid);
    free(writercmd);

    while (barcodes != NULL) {
        struct BarcodeInfo *bk;

        free(barcodes->name);
        free(barcodes->index);
        free(barcodes->delimiter);
        bk = barcodes->next;
        free(barcodes);
        barcodes = bk;
    }

    return r;
}
