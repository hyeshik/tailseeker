/*
 * sqi2fq.c
 *
 * Copyright (c) 2013 Hyeshik Chang
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

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#define SQILINEMAX      16384

static int
process(int cycle_begin, int cycle_end)
{
    char buf[SQILINEMAX];
    int use_istart=(cycle_begin == -1);

    cycle_begin--; /* to 0-based coordinate */

    for (;;) {
        char seqname[80];
        char *pch, *tile, *cluster, *seq, *qual, *istart;
        int i;

        if (fgets(buf, SQILINEMAX-1, stdin) == NULL)
            break;

        buf[SQILINEMAX-1] = 0;
        tile = cluster = seq = qual = istart = NULL;

        for (i = 0, pch = strtok(buf, "\t\r\n");
                pch != NULL; i++, pch = strtok(NULL, "\t\r\n"))
            switch (i) {
            case 0: /* tile */
                tile = pch;
                break;
            case 1: /* cluster */
                cluster = pch;
                break;
            case 2: /* qcpass */
                break;
            case 3: /* seq */
                seq = pch;
                break;
            case 4: /* qual */
                qual = pch;
                break;
            case 5: /* intensity */
                break;
            case 6: /* istart */
                istart = pch;
                break;
            default:
                break;
            }

        if (tile == NULL || cluster == NULL || seq == NULL || qual == NULL) {
            fprintf(stderr, "Invalid line found.\n");
            return 5;
        }

        snprintf(seqname, 79, "%s:%08d", tile, atoi(cluster));
        if (!use_istart) {
            seq += cycle_begin;
            qual += cycle_begin;
            seq[cycle_end - cycle_begin] = 0;
            qual[cycle_end - cycle_begin] = 0;
        }
        else {
            int r2start=atoi(istart);
            seq += r2start;
            qual += r2start;
        }

        printf("@%s\n%s\n+%s\n%s\n", seqname, seq, seqname, qual);
    }

    return 0;
}



int
main(int argc, char *argv[])
{
    int cycle_begin, cycle_end;

    if (argc < 3) {
        fprintf(stderr, "Usage: %s start end  (coordinates in 1-based inclusive)\n",
                        argv[0]);
        return 1;
    }

    cycle_begin = atoi(argv[1]);
    cycle_end = atoi(argv[2]);
    if (cycle_begin == -1 && cycle_end == -1)
        /* istart mode for read 2 */;
    else if (cycle_begin < 1 || cycle_end < 1 || cycle_begin > cycle_end) {
        fprintf(stderr, "Cycle range is wrong.\n");
        return 2;
    }

    return process(cycle_begin, cycle_end);
}
