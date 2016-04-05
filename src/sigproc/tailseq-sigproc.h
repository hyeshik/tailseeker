/*
 * tailseq-sigproc.h
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

#ifndef _TAILSEQ_SIGPROC_H_
#define _TAILSEQ_SIGPROC_H_

#include <stdio.h>
#include <stdint.h>
#include <zlib.h>


#define NUM_CHANNELS        4

#define PAFLAG_BARCODE_HAS_MISMATCHES       1
#define PAFLAG_DELIMITER_NOT_FOUND          2
#define PAFLAG_DELIMITER_HAS_MISMATCH       4
#define PAFLAG_DELIMITER_IS_SHIFTED         8
#define PAFLAG_BALANCER_BIASED              16
#define PAFLAG_BALANCER_SIGNAL_BAD          32
#define PAFLAG_MEASURED_FROM_FLUORESCENCE   64
#define PAFLAG_DARKCYCLE_EXISTS             128
#define PAFLAG_DARKCYCLE_OVER_THRESHOLD     256
#define PAFLAG_HAVE_3P_MODIFICATION         512


#define BCLREADER_OVERRIDDEN    ((struct BCLReader *)1)
                                /* placeholder for overridden cycles by an alternative call */
struct BCLReader {
    uint32_t nclusters;

    FILE *fptr;
    uint32_t read;
};

struct BCLData {
    size_t nclusters;
    uint8_t basequality[1];
};

struct IntensitySet {
    int16_t value[NUM_CHANNELS];
};

struct CIFReader {
    char filemagic[3];
    uint8_t version;
    uint8_t datasize;
    uint16_t first_cycle;
    uint16_t ncycles;
    uint32_t nclusters;

    FILE *fptr;
    uint32_t read;
    int16_t *readbuf;
    size_t readbuf_size;
};

struct CIFData {
    size_t nclusters;
    struct IntensitySet intensity[1];
};

struct UMIInterval;
struct UMIInterval {
    int start;
    int end;
    int length;
    int min_quality;
    float min_fraction_passes;
};

struct SampleInfo;
struct SampleInfo {
    char *name;
    char *index;

    char *delimiter;
    int delimiter_pos;
    int delimiter_length;

    char *fingerprint;
    int fingerprint_pos;
    int fingerprint_length;

    int maximum_index_mismatches;
    int maximum_delimiter_mismatches;
    int maximum_fingerprint_mismatches;

    struct UMIInterval *umi_ranges;
    int umi_ranges_count;
    int umi_total_length;

    FILE *stream;
    uint32_t clusters_mm0;
    uint32_t clusters_mm1;
    uint32_t clusters_mm2plus;
    uint32_t clusters_nodelim;
    uint32_t clusters_fpmismatch;
    uint32_t clusters_qcfailed;
    struct SampleInfo *next;
};

struct AlternativeCallReader {
    gzFile fptr;
    char *filename;
    int16_t ncycles;
    uint32_t read;
};

struct AlternativeCallInfo;
struct AlternativeCallInfo {
    char *filename;
    int first_cycle; /* the last cycle number will be determined from the file content */

    struct AlternativeCallReader *reader;
    struct AlternativeCallInfo *next;
};

#define CONTROL_NAME_MAX    128
struct ControlFilterInfo {
    char name[CONTROL_NAME_MAX];
    int first_cycle;
    int read_length;
    struct SampleInfo *barcode;
};

#define CONTROL_ALIGN_BASE_COUNT            5
#define CONTROL_ALIGN_MATCH_SCORE           1
#define CONTROL_ALIGN_MISMATCH_SCORE        2
#define CONTROL_ALIGN_GAP_OPEN_SCORE        4
#define CONTROL_ALIGN_GAP_EXTENSION_SCORE   1
#define CONTROL_ALIGN_MINIMUM_SCORE         0.65

struct PolyAFinderParameters {
    int weights[256];
    size_t max_terminal_modifications;
    size_t min_polya_length;
    size_t sigproc_trigger_polya_length;
};

#define T_INTENSITY_SCORE_BINS              200
struct PolyARulerParameters {
    int balancer_start;
    int balancer_end;
    int balancer_length;
    int balancer_minimum_occurrence;
    int balancer_num_positive_samples;
    int balancer_num_negative_samples;

    float colormatrix[NUM_CHANNELS * NUM_CHANNELS];

    float dark_cycles_threshold;
    int max_dark_cycles;

    float maximum_entropy;
    float t_intensity_score[T_INTENSITY_SCORE_BINS + 1];

    float polya_score_threshold;
    float downhill_extension_weight;
};

struct TailseekerConfig {
    /* section source */
    char *datadir, *laneid;
    int lane, tile, total_cycles;
    int index_start, index_length;
    int threep_start, threep_length;
    char *threep_colormatrix_filename;

    /* section options */
    int keep_no_delimiter;
    int threads;
    size_t read_buffer_size;
    int read_buffer_entry_count;

    /* section output */
    char *fastq_output;
    char *stats_output;
    char *length_dists_output;

    /* section alternative_calls */
    struct AlternativeCallInfo *altcalls;

    /* section control */
    struct ControlFilterInfo controlinfo;

    /* section polyA_finder */
    struct PolyAFinderParameters finderparams;

    /* section polyA_ruler */
    struct PolyARulerParameters rulerparams;
    float t_intensity_k, t_intensity_center;

    /* section sample: */
    struct SampleInfo *samples;
};


/* bclreader.c */
extern struct BCLReader *open_bcl_file(const char *filename);
extern void close_bcl_file(struct BCLReader *bcl);
extern int load_bcl_data(struct BCLReader *bcl, struct BCLData *data, uint32_t nclusters);
extern struct BCLData *new_bcl_data(uint32_t size);
extern void free_bcl_data(struct BCLData *data);
extern void format_basecalls(char *seq, char *qual, struct BCLData **basecalls,
                             int ncycles, uint32_t clusterno);
extern struct BCLReader **open_bcl_readers(const char *msgprefix, const char *datadir,
                                           int lane, int tile,
                                           int ncycles, struct AlternativeCallInfo *altcalls);
extern void close_bcl_readers(struct BCLReader **readers, int ncycles);

/* cifreader.c */
extern struct CIFReader *open_cif_file(const char *filename);
extern void close_cif_file(struct CIFReader *cif);
extern struct CIFData *new_cif_data(uint32_t size);
extern int load_cif_data(struct CIFReader *cif, struct CIFData *data, uint32_t nclusters);
extern void free_cif_data(struct CIFData *data);
extern void format_intensity(char *inten, struct CIFData **intensities,
                             int ncycles, uint32_t clusterno, int scalefactor);
extern void fetch_intensity(struct IntensitySet *signalout, struct CIFData **intensities,
                            int firstcycle, int ncycles, uint32_t clusterno);
extern struct CIFReader **open_cif_readers(const char *msgprefix, const char *datadir,
                                           int lane, int tile, int ncycles);
extern void close_cif_readers(struct CIFReader **readers, int ncycles);

/* altcalls.c */
extern struct AlternativeCallReader *open_alternative_calls(const char *filename);
extern int close_alternative_calls(struct AlternativeCallReader *acall, int checkend);
extern int load_alternative_calls(struct AlternativeCallReader *acall,
                                  struct BCLData **basecalls, uint32_t nclusters);
extern int open_alternative_calls_bundle(const char *msgprefix,
                                         struct AlternativeCallInfo *altcallinfo);
extern int close_alternative_calls_bundle(struct AlternativeCallInfo *altcallinfo,
                                          int checkend);

/* phix_control.c */
extern const char *phix_control_sequence;
extern const char *phix_control_sequence_rev;

/* controlaligner.c */
extern void initialize_ssw_score_matrix(int8_t *score_mat, int8_t match_score,
                                        int8_t mismatch_score);
extern ssize_t load_control_sequence(int8_t **control_seq);
extern int try_alignment_to_control(const char *sequence_read, const int8_t *control_seq,
                                    ssize_t control_seq_length,
                                    struct ControlFilterInfo *control_info,
                                    int8_t *ssw_score_mat, int32_t min_control_alignment_score,
                                    int32_t control_alignment_mask_len);

/* my_strstr.c */
extern char *my_strnstr(const char *s, const char *find, size_t len);

/* findpolya.c */
extern uint32_t find_polya(const char *seq, size_t seqlen,
                           struct PolyAFinderParameters *params);

/* signalproc.c */
extern int load_color_matrix(float *mtx, const char *filename);
extern int measure_polya_length(struct TailseekerConfig *cfg,
                struct CIFData **intensities, const char *sequence_formatted,
                uint32_t clusterno, int delimiter_end, int *procflags);
extern void precalc_score_tables(struct PolyARulerParameters *params,
                                 float t_score_k, float t_score_center);

/* spotanalyzer.c */
extern int process_spots(struct TailseekerConfig *cfg, uint32_t firstclusterno,
                         struct CIFData **intensities, struct BCLData **basecalls);

/* misc.c */
extern int inverse_4x4_matrix(const float *m, float *out);

/* parseconfig.c */
extern struct TailseekerConfig *parse_config(const char *filename);
extern void free_config(struct TailseekerConfig *cfg);

#endif
