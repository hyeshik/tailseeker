/*
 * utils.c
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
#include <stdarg.h>
#include "htslib/bgzf.h"


char *
replace_placeholder(const char *format, const char *old, const char *new)
{
    size_t newsize;
    char *found, *formatted, *wptr;

    found = strstr(format, old);
    if (found == NULL)
        return strdup(old);

    newsize = strlen(format) - strlen(old) + strlen(new);
    formatted = wptr = malloc(newsize + 1);
    strncpy(formatted, format, found - format);
    wptr += (found - format);
    strcpy(wptr, new);
    wptr += strlen(new);
    strcpy(wptr, found + strlen(old));

    return formatted;
}


int
bgzf_printf(BGZF *fp, const char *format, ...)
{
    va_list args;
    char buf[BUFSIZ];
    int r, len;

    va_start(args, format);
    len = vsprintf(buf, format, args);

    if (len > 0)
        r = bgzf_write(fp, buf, len);
    else if (len < 0)
        r = -1;
    else
        r = 0;

    va_end(args);

    return r;
}
