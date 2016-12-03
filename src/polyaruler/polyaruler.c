/*
 * polyaruler.c
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
#include <zlib.h>
#include <math.h>
#include <assert.h>
#include "../sigproc-flags.h"
#include "../signal-packs.h"


#define SCORE_LINEBUF_SIZE      16384
#define MAX_NUM_CYCLES          1024
#define TAGINFO_LINEBUF_SIZE    2048


static unpacked_score_t *
load_score_cutoffs(const char *filename, const char *tileid,
                   ssize_t *ncycles)
{
    FILE *fp;
    char linebuf[SCORE_LINEBUF_SIZE], *bptr, *cutofftok;
    size_t tileid_len;
    unpacked_score_t *cutoffs;
    int i;

    tileid_len = strlen(tileid);
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s.\n", filename);
        return NULL;
    }

    /* Search for the line for the designated tile. */
    for (;;) {
        if (fgets(linebuf, sizeof(linebuf), fp) == NULL) {
            fprintf(stderr, "No poly(A) score cutoffs found for tile %s.\n",
                    tileid);
            fclose(fp);
            return NULL;
        }

        if (strlen(linebuf) >= tileid_len &&
                memcmp(linebuf, tileid, tileid_len) == 0 &&
                linebuf[tileid_len] == '\t')
            break;
    }

    fclose(fp);

    cutoffs = malloc(sizeof(unpacked_score_t) * MAX_NUM_CYCLES);
    if (cutoffs == NULL)
        return NULL;

    bptr = linebuf + tileid_len + 1;
    for (i = 0; (cutofftok = strsep(&bptr, "\t\r\n, ")) != NULL; i++) {
        if (*cutofftok == 0)
            break;

        if (i >= MAX_NUM_CYCLES) {
            fprintf(stderr, "Exceeded the maximum allowed number of cycles. "
                            "Adjust MAX_NUM_CYCLES in " __FILE__ ".\n");
            free(cutoffs);
            return NULL;
        }

        cutoffs[i] = 1 + (int)(atof(cutofftok) *
                               (double)(SIGNALPACKET_SCORE_MAX - 1));
        if (cutoffs[i] >= SIGNALPACKET_SCORE_MAX)
            cutoffs[i] = SIGNALPACKET_SCORE_MAX;
    }

    *ncycles = i;

    return cutoffs;
}

static int
measure_polya_length(const signal_packet_t *read_scores,
                     ssize_t read_num_cycles,
                     int first_cycle,
                     const unpacked_score_t *score_cutoffs,
                     ssize_t cutoffs_num_cycles,
                     float downhill_ext_weight)
{
    ssize_t i, physical_cycle;
    int score_cumsum, score_cumsum_max, score_argmax_cycle;
    int downhill_ext;

    score_cumsum = score_cumsum_max = 0;
    score_argmax_cycle = -1;
    physical_cycle = first_cycle;

    for (i = 0; i < read_num_cycles && physical_cycle < cutoffs_num_cycles;
            i++, physical_cycle++) {
        unpacked_score_t cutoff = score_cutoffs[physical_cycle];

        if (cutoff == 0) /* NaN in the integer representation */
            continue;

        score_cumsum += -1 + (read_scores[i].score >= cutoff) * 2;
        if (score_cumsum > score_cumsum_max) {
            score_cumsum_max = score_cumsum;
            score_argmax_cycle = i;
        }
    }

    downhill_ext = -1;

    /* Try extending the poly(A) until entropy continuously decreases. */
    if (score_argmax_cycle >= 0)
        for (i = score_argmax_cycle + 1; i < read_num_cycles; i++)
            if (read_scores[i].score > 0) { /* Non-dark signals */
                if (read_scores[i].downhill > 0)
                    downhill_ext = i;
                else
                    break;
            }

    if (downhill_ext >= 0)
        return score_argmax_cycle + 1 +
            (int)roundf((downhill_ext - score_argmax_cycle) *
                        downhill_ext_weight);
    else
        return score_argmax_cycle + 1;
}

static void
add_polya_score_sample(cluster_count_t *samplecounts,
                       signal_packet_t *scores, int length,
                       int firstcycle, int sampling_bins)
{
    cluster_count_t *sptr;
    int i, binno;

    sptr = samplecounts + (firstcycle * sampling_bins);
    for (i = 0; i < length; i++, scores++, sptr += sampling_bins)
        if (scores->score > 0) {
            binno = ((scores->score - 1.) / (SIGNALPACKET_SCORE_MAX - 1)) *
                    sampling_bins;
            if (binno >= sampling_bins)
                binno = sampling_bins - 1;

            sptr[binno]++;
        }   
}

static int 
write_signal_samples_dists(const char *filename,
                           const cluster_count_t *counts,
                           int total_cycles, int sampling_bins)
{
    uint32_t header_elements[3];
    gzFile fp;
    int r;

    fp = gzopen(filename, "wb");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s to write.\n", filename);
        return -1; 
    }   

    header_elements[0] = sizeof(cluster_count_t);
    header_elements[1] = total_cycles;
    header_elements[2] = sampling_bins;

    r = (gzwrite(fp, header_elements, sizeof(header_elements)) < 0 ||
         gzwrite(fp, counts, sizeof(cluster_count_t) *
                    total_cycles * sampling_bins) < 0) ? -1 : 0;
    gzclose(fp);

    return r;
}

static int16_t *
process_polya_ruling(const char *filename,
                     const unpacked_score_t *score_cutoffs,
                     ssize_t cutoffs_num_cycles, int minimum_polya_len,
                     float downhill_ext_weight, int dist_sampling_bins,
                     float dist_sampling_gap, const char *sigdist_output,
                     size_t *ret_elements)
{
    gzFile fp;
    ssize_t record_size, bytesread, sigdist_size;
    int16_t *polya_measurements;
    int r;
    cluster_count_t *pos_score_counts;
    struct {
        uint32_t elemsize;
        uint32_t total_clusters;
        uint32_t max_cycles;
    } header;

    sigdist_size = cutoffs_num_cycles * dist_sampling_bins * sizeof(cluster_count_t);
    pos_score_counts = malloc(sigdist_size);
    if (pos_score_counts == NULL) {
        perror("process_polya_ruling");
        return NULL;
    }
    memset(pos_score_counts, 0, sigdist_size);

    fp = gzopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open the input file: %s\n", filename);
        free(pos_score_counts);
        return NULL;
    }

    if (gzread(fp, &header, sizeof(header)) != sizeof(header)) {
        fprintf(stderr, "Failed to read the header from %s.\n", filename);
        gzclose(fp);
        free(pos_score_counts);
        return NULL;
    }

    if (header.elemsize != sizeof(signal_packet_t)) {
        fprintf(stderr, "The file %s was written in a machine with "
                        "different architecture.\n", filename);
        gzclose(fp);
        free(pos_score_counts);
        return NULL;
    }

    polya_measurements = malloc(sizeof(int16_t) * header.total_clusters);
    if (polya_measurements == NULL) {
        perror("process_polya_ruling");
        gzclose(fp);
        free(pos_score_counts);
        return NULL;
    }

    memset(polya_measurements, 0xff,
           sizeof(int16_t) * header.total_clusters); /* fill with -1 */
    record_size = sizeof(struct SignalRecordHeader) +
                  header.max_cycles * sizeof(signal_packet_t);

    for (;;) {
        int polya_len;
        struct TSRecord {
            struct SignalRecordHeader header;
            signal_packet_t scores[];
        } rec;

        bytesread = gzread(fp, (void *)&rec, record_size);
        if (bytesread == 0)
            break;

        if (bytesread < record_size) {
            fprintf(stderr, "Unexpected end of file.\n");
            gzclose(fp);
            free(pos_score_counts);
            return NULL;
        }

        polya_len = measure_polya_length(rec.scores,
                rec.header.valid_cycle_count, rec.header.first_cycle,
                score_cutoffs, cutoffs_num_cycles, downhill_ext_weight);

        if (polya_len >= minimum_polya_len) {
            int sampling_len;
            polya_measurements[rec.header.clusterno] = polya_len;
            sampling_len = polya_len - (int)(polya_len * dist_sampling_gap);
            add_polya_score_sample(pos_score_counts, rec.scores,
                sampling_len, rec.header.first_cycle, dist_sampling_bins);
        }
    }

    gzclose(fp);

    r = write_signal_samples_dists(sigdist_output, pos_score_counts,
                    cutoffs_num_cycles, dist_sampling_bins);
    free(pos_score_counts);
    if (r < 0)
        return NULL;

    *ret_elements = header.total_clusters;

    return polya_measurements;
}

static int
output_corrected_polya_measurements(const char *taginfo_file,
                                    int16_t *polya_measurements,
                                    size_t nclusters,
                                    const char *tile_id)
{
    gzFile fp;

    fp = gzopen(taginfo_file, "rt");
    if (fp == NULL)
        return -1;

    for (;;) {
        char linebuf[TAGINFO_LINEBUF_SIZE];
        char *bptr, *token, *mods, *umi;
        uint32_t clusterno, flags;
        int i;

        if (gzgets(fp, linebuf, TAGINFO_LINEBUF_SIZE) == NULL) {
            int errno;
            (void)gzerror(fp, &errno);
            if (errno == 0)
                break; /* end-of-file */
            else {
                fprintf(stderr, "Error occurred on reading %s.\n",
                        taginfo_file);
                gzclose(fp);
                return -1;
            }
        }

        bptr = linebuf;
        token = strsep(&bptr, "\t\n\r");
        if (token == NULL)
            continue;

        clusterno = atoi(token);
        if (polya_measurements[clusterno] < 0) {
            /* Poly(A) length is not revised. Bypass the line. */
            token[strlen(token)] = '\t';
            fputs(tile_id, stdout);
            fputc('\t', stdout);
            fputs(linebuf, stdout);
            continue;
        }

        mods = umi = NULL;
        flags = 0;

        for (i = 0; (token = strsep(&bptr, "\t\n\r")) != NULL; i++)
            switch (i) {
            case 0: flags = atoi(token); break;
            case 1: /*polya_prelim = atoi(token);*/ break;
            case 2: mods = token; break;
            case 3: umi = token; break;
            default: break;
            }

        flags |= PAFLAG_MEASURED_FROM_FLUORESCENCE;
        printf("%s\t%d\t%d\t%d\t%s\t%s\n", tile_id, clusterno, flags,
               polya_measurements[clusterno], mods, umi);
    }

    return 0;
}

int
main(int argc, char *argv[])
{
    ssize_t cutoffs_num_cycles;
    unpacked_score_t *score_cutoffs;
    int16_t *polya_measurements;
    size_t nclusters;
    const char *tile_id, *signals_file, *cutoffs_file, *taginfo_file;
    const char *sigdist_file;
    float downhill_ext_weight, dist_sampling_gap;
    int minimum_polya_len, dist_sampling_bins;

    if (argc < 10) {
        fprintf(stderr, "Usage: %s {tile id} {signals} {cutoffs} {min polya} "
                        "{downhill ext weight} {taginfo} {sampling bin count} "
                        "{positive sampling gap} {positive sampling output}\n",
                        argv[0]);
        return 1;
    }

    tile_id = argv[1];
    signals_file = argv[2];
    cutoffs_file = argv[3];
    minimum_polya_len = atoi(argv[4]);
    downhill_ext_weight = atof(argv[5]);
    taginfo_file = argv[6];
    dist_sampling_bins = atoi(argv[7]);
    dist_sampling_gap = atof(argv[8]);
    sigdist_file = argv[9];

    /* Load per-cycle poly(A) score cutoffs table */
    score_cutoffs = load_score_cutoffs(cutoffs_file, tile_id, &cutoffs_num_cycles);
    if (score_cutoffs == NULL)
        return 1;

    /* Measure the lengths of poly(A) tails */
    polya_measurements = process_polya_ruling(signals_file, score_cutoffs,
            cutoffs_num_cycles, minimum_polya_len, downhill_ext_weight,
            dist_sampling_bins, dist_sampling_gap, sigdist_file,
            &nclusters);
    free(score_cutoffs);
    if (polya_measurements == NULL)
        return 2;

    /* Apply the measurements to the existing taginfo */
    if (output_corrected_polya_measurements(taginfo_file,
                polya_measurements, nclusters, tile_id) < 0) {
        free(polya_measurements);
        return 3;
    }

    free(polya_measurements);

    return 0;
}
