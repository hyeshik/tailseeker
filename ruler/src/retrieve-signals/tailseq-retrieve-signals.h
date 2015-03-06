/*
 * tailseq-retrieve-signals.h
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

#ifndef _TAILSEQ_RETRIEVE_SIGNALS_H_
#define _TAILSEQ_RETRIEVE_SIGNALS_H_

#define NUM_CHANNELS    4

struct BCLHandler {
    uint16_t ncycles;
    uint32_t nclusters;

    FILE *fptr;
    uint32_t read;
};

struct BCLData {
    size_t nclusters;
    uint8_t *base;
    uint8_t *quality;
};

struct IntensitySet {
    uint16_t value[NUM_CHANNELS];
};

struct CIFHandler {
    char filemagic[3];
    uint8_t version;
    uint8_t datasize;
    uint16_t first_cycle;
    uint16_t ncycles;
    uint32_t nclusters;

    FILE *fptr;
    uint32_t read;
};

struct CIFData {
    size_t nclusters;
    struct IntensitySet intensity[1];
};



/* bclreader.c */
extern struct BCLData *load_bcl_file(const char *filename);
extern void free_bcl_data(struct BCLData *data);
extern void format_basecalls(char *seq, char *qual, struct BCLData **basecalls,
                             int ncycles, uint32_t clusterno);


/* cifreader.c */
extern struct CIFData *load_cif_file(const char *filename);
extern void free_cif_data(struct CIFData *data);
extern void format_intensity(char *inten, struct CIFData **intensities,
                             int ncycles, uint32_t clusterno, double scalefactor);

/* phix_control.c */
extern const char *phix_control_sequence;

#endif
