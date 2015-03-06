/*
 * cifreader.c
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
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <endian.h>
#include "tailseq-retrieve-signals.h"


#define BATCH_BLOCK_SIZE    65536

static const char BASE64_ENCODE_TABLE[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef"
                                          "ghijklmnopqrstuvwxyz0123456789+/";


struct CIFData *
load_cif_file(const char *filename)
{
    FILE *fp;
    struct CIFHandler header;
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


void
free_cif_data(struct CIFData *data)
{
    free(data);
}


void
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
