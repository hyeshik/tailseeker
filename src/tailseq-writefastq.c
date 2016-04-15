/*
 * tailseq-writefastq.c
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
#include <getopt.h>
#include <zlib.h>
#include <htslib/bgzf.h>
#include "utils.h"

#define GZIP_READ_BUFFER_SIZE   1024*1024

#define LINE_BUFFER_SIZE        8192
#define MAX_TILENAME_LEN        63
#define MAX_MODIFICATION_LEN    50
#define MAX_ENTRYNAME_LEN       255

struct TagInfo {
    char tilename[MAX_TILENAME_LEN+1];
    int clusterno;
    int flags;
    int polyA_len;

    char modifications[MAX_MODIFICATION_LEN+1];
    size_t modification_len;

    int num_duplicates;
};


static inline int
parse_taginfo_line(struct TagInfo *taginfo, const char *line)
{
    size_t length;
    char *pos, *start;

    /* Locate tilename and copy */
    pos = strchr(line, '\t');
    if (pos == NULL)
        return -1;
    length = (size_t)(pos - line);
    memcpy(taginfo->tilename, line, length);
    taginfo->tilename[length] = '\0';

    /* Locate cluster number and parse it */
    start = pos + 1;
    taginfo->clusterno = strtoul(start, &pos, 10);
    if (pos == NULL)
        return -1;

    /* Locate flags and parse it */
    start = pos + 1;
    taginfo->flags = strtol(start, &pos, 10);
    if (pos == NULL)
        return -1;

    /* Locate poly(A) length and parse it */
    start = pos + 1;
    taginfo->polyA_len = strtol(start, &pos, 10);
    if (pos == NULL)
        return -1;

    /* Locate 3' end modification sequence and copy */
    start = pos + 1;
    pos = strchr(start, '\t');
    if (pos == NULL)
        return -1;
    length = (size_t)(pos - start);
    if (length > 0)
        memcpy(taginfo->modifications, start, length);
    taginfo->modifications[length] = '\0';
    taginfo->modification_len = length;

    /* Locate duplicity and parse it */
    start = pos + 1;
    taginfo->num_duplicates = strtol(start, &pos, 10);
    if (pos == NULL)
        return -1;

    return 0;
}

static int
exchange_seqqual_handle(gzFile *fp, const char *filename_format, const char *tileid)
{
    char *filename;

    if (*fp != NULL)
        gzclose(*fp);

    filename = replace_placeholder(filename_format, "@tile@", tileid);
    if (filename == NULL) {
        perror("exchange_seqqual_handle");
        return -1;
    }

    *fp = gzopen(filename, "rt");
    if (*fp == NULL) {
        fprintf(stderr, "Failed to open %s.\n", filename);
        free(filename);
        return -1;
    }

    free(filename);
    return 0;
}

static inline int
read_and_peek_seqqual(gzFile fp, char *line)
{
    char *endptr;
    int clusterno;

    if (gzgets(fp, line, LINE_BUFFER_SIZE) == NULL)
        return -1;

    clusterno = (int)strtol(line, &endptr, 10);
    if (endptr == NULL || *endptr != '\t')
        return -1;

    return clusterno;
}

static int
write_fastq_entry(BGZF *fastq5out, BGZF *fastq3out,
                  struct TagInfo *taginfo, const char *seqqual_line)
{
    char entryname[MAX_ENTRYNAME_LEN];
    char buf[LINE_BUFFER_SIZE], *bufp, *sqptr_1, *sqptr_2;
    int entryname_len, seq_len;

    entryname_len = snprintf(entryname, MAX_ENTRYNAME_LEN, "%s:%08u:%04x:%d:%d:%s",
                           taginfo->tilename, (unsigned)taginfo->clusterno,
                           taginfo->flags, taginfo->num_duplicates, taginfo->polyA_len,
                           taginfo->modifications);

    bufp = buf;

    /* Line 1 of the 5'-side FASTQ */
    *bufp++ = '@';
    memcpy(bufp, entryname, entryname_len);
    bufp += entryname_len;
    *bufp++ = '\n';

    /* Line 2 of the 5'-side FASTQ */
    /* - find the tab between cluster no. and 5' seq */
    sqptr_1 = strchr(seqqual_line, '\t');
    if (sqptr_1 == NULL)
        return -1;

    /* - find the tab between 5' seq and 5' qual */
    sqptr_2 = strchr(++sqptr_1, '\t');
    if (sqptr_2 == NULL)
        return -1;

    seq_len = (int)(sqptr_2 - sqptr_1);
    memcpy(bufp, sqptr_1, seq_len);
    bufp += seq_len;
    *bufp++ = '\n';

    /* Line 3 of the 5'-side FASTQ */
    *bufp++ = '+';
    memcpy(bufp, entryname, entryname_len);
    bufp += entryname_len;
    *bufp++ = '\n';

    /* Line 4 of the 5'-side FASTQ */
    sqptr_1 = strchr(++sqptr_2, '\t');
    if (sqptr_1 == NULL)
        return -1;
    if (seq_len != (int)(sqptr_1 - sqptr_2))
        return -1;
    memcpy(bufp, sqptr_2, seq_len);
    bufp += seq_len;
    *bufp++ = '\n';

    if (bgzf_write(fastq5out, buf, (size_t)(bufp - buf)) < 0)
        return -1;


    bufp = buf;

    /* Line 1 of the 3'-side FASTQ */
    *bufp++ = '@';
    memcpy(bufp, entryname, entryname_len);
    bufp += entryname_len;
    *bufp++ = '\n';

    /* Line 2 of the 3'-side FASTQ */
    sqptr_2 = strchr(++sqptr_1, '\t');
    if (sqptr_2 == NULL)
        return -1;
    seq_len = (int)(sqptr_2 - sqptr_1);
    memcpy(bufp, sqptr_1, seq_len);
    bufp += seq_len;
    *bufp++ = '\n';

    /* Line 3 of the 3'-side FASTQ */
    *bufp++ = '+';
    memcpy(bufp, entryname, entryname_len);
    bufp += entryname_len;
    *bufp++ = '\n';

    /* Line 4 of the 3'-side FASTQ */
    sqptr_1 = strchr(++sqptr_2, '\n');
    if (sqptr_1 == NULL)
        return -1;
    if (seq_len != (int)(sqptr_1 - sqptr_2))
        return -1;
    memcpy(bufp, sqptr_2, seq_len);
    bufp += seq_len;
    *bufp++ = '\n';

    return bgzf_write(fastq3out, buf, (size_t)(bufp - buf));
}

static int
process_write_fastq(gzFile taginfof, const char *seqqual_filename,
                    BGZF *fastq5out, BGZF *fastq3out)
{
    char tile_current[MAX_TILENAME_LEN];
    char line[LINE_BUFFER_SIZE];
    struct TagInfo taginfo;
    int errnum;
    const char *message;
    gzFile seqqualf;

    tile_current[0] = '\0';
    seqqualf = NULL;

    for (; gzgets(taginfof, line, LINE_BUFFER_SIZE) != NULL;) {
        int seqqual_clusterno;

        if (parse_taginfo_line(&taginfo, line) < 0) {
            fprintf(stderr, "Failed to parse a line: %s", line);
            return -1;
        }

        /* On tile location changes, reopen the seqqual file. */
        if (tile_current[0] == '\0' || strcmp(tile_current, taginfo.tilename) != 0) {
            strcpy(tile_current, taginfo.tilename);

            if (exchange_seqqual_handle(&seqqualf, seqqual_filename, tile_current) < 0)
                return -1;
        }

        /* Read lines from the seqqual file until cluster number matches to taginfo's. */
        do {
            seqqual_clusterno = read_and_peek_seqqual(seqqualf, line);

            if (seqqual_clusterno < 0) {
                fprintf(stderr, "Unexpected format in a seqqual file.\n");
                gzclose(seqqualf);
                return -1;
            }
            else if (seqqual_clusterno > taginfo.clusterno) {
                fprintf(stderr, "The taginfo file must be a subset of seqqual.\n");
                gzclose(seqqualf);
                return -1;
            }
        } while (seqqual_clusterno < taginfo.clusterno);

        if (write_fastq_entry(fastq5out, fastq3out, &taginfo, line) < 0) {
            fprintf(stderr, "Failed to write an entry in %s.\n", tile_current);
            gzclose(seqqualf);
            return -1;
        }
    }

    message = gzerror(taginfof, &errnum);
    if (errnum != 0) {
        fprintf(stderr, "process_write_fastq: %s\n", message);
        return -1;
    }

    if (seqqualf != NULL)
        gzclose(seqqualf);

    return 0;
}


int
main(int argc, char *argv[])
{
    gzFile taginfof;
    BGZF *fastq5, *fastq3;
    int r, threads;
    char *taginfo_filename, *seqqual_filename;
    char *fastq5_filename, *fastq3_filename;

    struct option long_options[] =
    {
        {"threads", required_argument,  0,  't'},
        {"taginfo", required_argument,  0,  'i'},
        {"seqqual", required_argument,  0,  's'},
        {"fastq5",  required_argument,  0,  '5'},
        {"fastq3",  required_argument,  0,  '3'},
        {0, 0, 0, 0}
    };

    threads = 1;
    taginfo_filename = seqqual_filename =
        fastq5_filename = fastq3_filename = NULL;

    while (1) {
        int option_index=0;
        int c;

        c = getopt_long(argc, argv, "t:i:s:5:3:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
            case 't': /* --threads */
                threads = atoi(optarg);
                break;

            case 'i': /* --taginfo */
                taginfo_filename = strdup(optarg);
                break;

            case 's': /* --seqqual */
                seqqual_filename = strdup(optarg);
                break;

            case '5': /* --fastq5 */
                fastq5_filename = strdup(optarg);
                break;

            case '3': /* --fastq3 */
                fastq3_filename = strdup(optarg);
                break;

            default:
                abort();
        }
    }

    if (taginfo_filename == NULL || seqqual_filename == NULL ||
        fastq5_filename == NULL || fastq3_filename == NULL) {
        fprintf(stderr, "One or more required arguments are not specified.\n");
        return -1;
    }

    taginfof = gzopen(taginfo_filename, "rt");
    if (taginfof == NULL) {
        fprintf(stderr, "Failed to open %s.\n", taginfo_filename);
        return -1;
    }

    gzbuffer(taginfof, GZIP_READ_BUFFER_SIZE);

    fastq5 = bgzf_open(fastq5_filename, "w");
    if (fastq5 == NULL) {
        fprintf(stderr, "Failed to open %s.\n", fastq5_filename);
        gzclose(taginfof);
        return -1;
    }
    if (threads > 1 && bgzf_mt(fastq5, threads, 256) < 0) {
        fprintf(stderr, "Failed to setting multi-threading for FASTQ5\b");
        return -1;
    }

    fastq3 = bgzf_open(fastq3_filename, "w");
    if (fastq3 == NULL) {
        fprintf(stderr, "Failed to open %s.\n", fastq5_filename);
        gzclose(taginfof);
        bgzf_close(fastq5);
        return -1;
    }
    if (threads > 1 && bgzf_mt(fastq3, threads, 256) < 0) {
        fprintf(stderr, "Failed to setting multi-threading for FASTQ3\b");
        return -1;
    }

    r = process_write_fastq(taginfof, seqqual_filename, fastq5, fastq3);

    gzclose(taginfof);
    bgzf_close(fastq5);
    bgzf_close(fastq3);

    return r;
}
