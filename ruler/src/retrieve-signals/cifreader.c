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


static const char BASE64_ENCODE_TABLE[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef"
                                          "ghijklmnopqrstuvwxyz0123456789+/";


struct CIFReader *
open_cif_file(const char *filename)
{
    FILE *fp;
    struct CIFReader *hdl;

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "open_cif_file: Can't open file %s.\n", filename);
        return NULL;
    }

    hdl = malloc(sizeof(struct CIFReader));
    if (hdl == NULL) {
        perror("open_cif_file");
        return NULL;
    }

    if (fread(hdl->filemagic, 3, 1, fp) < 1)
        goto onError;

    if (memcmp(hdl->filemagic, "CIF", 3) != 0) {
        fprintf(stderr, "File %s does not seem to be a CIF file.", filename);
        goto onError;
    }

    /* read version and data size */
    if (fread(&hdl->version, 2, 1, fp) < 1) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        goto onError;
    }
    if (hdl->version != 1) {
        fprintf(stderr, "Unsupported CIF version: %d.\n", hdl->version);
        goto onError;
    }
    if (hdl->datasize != 2) {
        fprintf(stderr, "Unsupported data size (%d).\n", hdl->datasize);
        goto onError;
    }

    /* read first cycle and number of cycles */
    if (fread(&hdl->first_cycle, 2, 2, fp) < 2) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        goto onError;
    }
    if (fread(&hdl->nclusters, 4, 1, fp) < 1) {
        fprintf(stderr, "Unexpected EOF %s:%d.\n", __FILE__, __LINE__);
        goto onError;
    }

    hdl->first_cycle = le16toh(hdl->first_cycle);
    hdl->ncycles = le16toh(hdl->ncycles);
    hdl->nclusters = le32toh(hdl->nclusters);
    hdl->fptr = fp;
    hdl->read = 0;

    return hdl;

  onError:
    free(hdl);
    fclose(fp);
    return NULL;
}


void
close_cif_file(struct CIFReader *cif)
{
    if (cif->fptr != NULL)
        fclose(cif->fptr);
    free(cif);
}


int
load_cif_data(struct CIFReader *cif, struct CIFData *data, uint32_t nclusters)
{
    uint32_t toread;

    if (cif->read >= cif->nclusters) {
        data->nclusters = 0;
        return 0;
    }

    if (cif->read + nclusters >= cif->nclusters) /* does file has enough clusters to read? */
        toread = cif->nclusters - cif->read; /* all clusters left */
    else
        toread = nclusters;

    data->nclusters = toread;
    if (fread(data->intensity, sizeof(struct IntensitySet), toread, cif->fptr) != toread) {
        fprintf(stderr, "Not all data were loaded from a CIF.");
        return -1;
    }

    cif->read += toread;

    return 0;
}


struct CIFData *
new_cif_data(uint32_t size)
{
    struct CIFData *cdata;

    cdata = malloc(sizeof(size_t) + sizeof(struct IntensitySet) * size);
    if (cdata == NULL) {
        perror("new_cif_data");
        return NULL;
    }

    cdata->nclusters = 0;
    return cdata;
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


struct CIFReader **
open_cif_readers(const char *msgprefix, const char *datadir, int lane, int tile, int ncycles)
{
    int cycleno;
    char path[PATH_MAX];
    struct CIFReader **readers;

    readers = malloc(sizeof(struct CIFReader *) * ncycles);
    if (readers == NULL) {
        perror("open_cif_readers");
        return NULL;
    }

    printf("%sUsing cluster intensities from %s.\n", msgprefix, datadir);

    for (cycleno = 0; cycleno < ncycles; cycleno++) {
        snprintf(path, PATH_MAX, "%s/L%03d/C%d.1/s_%d_%04d.cif", datadir, lane, cycleno+1,
                 lane, tile);

        readers[cycleno] = open_cif_file(path);
        if (readers[cycleno] == NULL) {
            while (--cycleno >= 0)
                close_cif_file(readers[cycleno]);
            return NULL;
        }
    }

    return readers;
}

void
close_cif_readers(struct CIFReader **readers, int ncycles)
{
    int cycleno;

    for (cycleno = 0; cycleno < ncycles; cycleno++)
        close_cif_file(readers[cycleno]);
}
