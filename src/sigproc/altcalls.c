/*
 * altcalls.c
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
#include <errno.h>
#include <string.h>
#include <limits.h>
#include "tailseq-sigproc.h"


static const uint8_t DNABASE2NUM_ac[128] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};


static int16_t
check_fastq_read_length(gzFile hdl, const char *filename)
{
    char buf[BUFSIZ];

    if (gzgets(hdl, buf, BUFSIZ) == NULL) {
        fprintf(stderr, "No sequence available in %s.\n", filename);
        return -1;
    }

    if (buf[0] != '@') {
        fprintf(stderr, "%s is not in FASTQ format.\n", filename);
        return -1;
    }

    if (gzgets(hdl, buf, BUFSIZ) == NULL) {
        fprintf(stderr, "Unexpected EOF while reading %s\n", filename);
        return -1;
    }

    if (gzrewind(hdl) < 0) {
        fprintf(stderr, "Failed to rewind the FASTQ reader for %s\n", filename);
        return -1;
    }

    return strlen(buf) - 1;
}


struct AlternativeCallReader *
open_alternative_calls(const char *filename)
{
    struct AlternativeCallReader *acall;
    gzFile hdl;

    hdl = gzopen(filename, "rt");
    if (hdl == NULL) {
        perror("open_alternative_calls");
        return NULL;
    }

    acall = malloc(sizeof(struct AlternativeCallReader));
    if (acall == NULL) {
        perror("open_alternative_calls");
        gzclose(hdl);
        return NULL;
    }

    acall->fptr = hdl;
    acall->ncycles = check_fastq_read_length(hdl, filename);
    if (acall->ncycles < 0) {
        perror("open_alternative_calls");
        close_alternative_calls(acall, 0);
        return NULL;
    }

    acall->read = 0;
    acall->filename = strdup(filename);
    if (acall->filename == NULL) {
        perror("open_alternative_calls");
        close_alternative_calls(acall, 0);
        return NULL;
    }

    return acall;
}


int
close_alternative_calls(struct AlternativeCallReader *acall, int checkend)
{
    char buf[BUFSIZ];

    if (checkend && gzgets(acall->fptr, buf, BUFSIZ) != NULL) {
        fprintf(stderr, "Extra sequences found in %s\n", acall->filename);
        return -1;
    }

    if (acall->filename != NULL)
        free(acall->filename);

    if (acall->fptr != NULL)
        gzclose(acall->fptr);

    free(acall);

    return 0;
}


int
load_alternative_calls(struct AlternativeCallReader *acall, struct BCLData **basecalls,
                       uint32_t nclusters)
{
    char buf[BUFSIZ], seqbuf[BUFSIZ], qualbuf[BUFSIZ];
    uint32_t clusterno;
    int16_t ncycles, j;

    ncycles = acall->ncycles;

    for (clusterno = 0; clusterno < nclusters; clusterno++) {
        size_t seqlen;

        /* header 1 */
        if (gzgets(acall->fptr, buf, BUFSIZ) == NULL) {
            fprintf(stderr, "Not enough sequences from %s\n", acall->filename);
            return -1;
        }

        if (buf[0] != '@') {
            fprintf(stderr, "%s is not in FASTQ format.\n", acall->filename);
            return -1;
        }

        /* sequence */
        if (gzgets(acall->fptr, seqbuf, BUFSIZ) == NULL) {
            fprintf(stderr, "Unexpected EOF while reading %s\n", acall->filename);
            return -1;
        }

        seqlen = strlen(seqbuf) - 1;
        if (ncycles == -1)
            ncycles = acall->ncycles = seqlen;
        else if (seqlen != ncycles) {
            fprintf(stderr, "Sequence length in the FASTQ is inconsistent.\n");
            return -1;
        }

        /* header 2 */
        if (gzgets(acall->fptr, buf, BUFSIZ) == NULL) {
            fprintf(stderr, "Unexpected EOF while reading %s\n", acall->filename);
            return -1;
        }

        if (buf[0] != '+') {
            fprintf(stderr, "%s is not in FASTQ format.\n", acall->filename);
            return -1;
        }

        /* quality */
        if (gzgets(acall->fptr, qualbuf, BUFSIZ) == NULL) {
            fprintf(stderr, "Unexpected EOF while reading %s\n", acall->filename);
            return -1;
        }

        seqlen = strlen(qualbuf) - 1;
        if (seqlen != ncycles) {
            fprintf(stderr, "Sequence length in the FASTQ is inconsistent.\n");
            return -1;
        }

        for (j = 0; j < ncycles; j++)
            basecalls[j]->basequality[clusterno] = (
                DNABASE2NUM_ac[(int)seqbuf[j]] | (((uint8_t)(qualbuf[j] - 33)) << 2));
    }

    for (j = 0; j < ncycles; j++)
        basecalls[j]->nclusters = nclusters;

    return 0;
}


int
open_alternative_calls_bundle(const char *msgprefix, struct AlternativeCallInfo *altcallinfo)
{
    for (; altcallinfo != NULL; altcallinfo = altcallinfo->next)
        if (altcallinfo->reader == NULL) {
            altcallinfo->reader = open_alternative_calls(altcallinfo->filename);
            if (altcallinfo->reader == NULL)
                return -1;

            printf("%sUsing FASTQ %s for cycles %d-%d.\n", msgprefix,
                    altcallinfo->filename, altcallinfo->first_cycle + 1,
                    altcallinfo->first_cycle + altcallinfo->reader->ncycles);
        }

    return 0;
}


int
close_alternative_calls_bundle(struct AlternativeCallInfo *altcallinfo, int checkend)
{
    for (; altcallinfo != NULL; altcallinfo = altcallinfo->next)
        if (altcallinfo->reader != NULL &&
                close_alternative_calls(altcallinfo->reader, checkend) == -1)
            return -1;

    return 0;
}

