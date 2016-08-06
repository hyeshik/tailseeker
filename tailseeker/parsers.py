#!/usr/bin/env python3
#
# Copyright (c) 2016 Hyeshik Chang
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# - Hyeshik Chang <hyeshik@snu.ac.kr>
#

from .fileutils import LineParser, open_gzip_buffered

decode_bytes = lambda x: x.decode()

parse_taginfo = LineParser([
    ('tile', None),
    ('cluster', int),
    ('pflags', int),
    ('polyA', int),
    ('mods', None),
    ('clones', int)
], linefeed=b'\n')

parse_taginfo_internal = LineParser([
    ('tile', None),
    ('cluster', int),
    ('pflags', int),
    ('polyA', int),
    ('mods', None),
    ('umi', None)
], linefeed=b'\n')

parse_refined_taginfo = LineParser([
    ('tile', None),
    ('cluster', int),
    ('pflags', int),
    ('clones', int),
    ('polyA', int),
    ('unaligned_polyA', int),
    ('mods', None),
    ('unaligned_mods', None),
], linefeed=b'\n')

parse_sam = LineParser([
    ('qname', None),
    ('flag', int),
    ('rname', None),
    ('pos', int), # caution! A SAM file uses 1-based coordinate system.
    ('mapq', int),
    ('cigar', None),
    ('rnext', None),
    ('pnext', None),
    ('tlen', int),
    ('seq', None),
    ('qual', None),
], comment=b'@')

def parse_fastq(filename):
    linebuf = []
    for line in open_gzip_buffered(filename, 'rb'):
        linebuf.append(line)
        if len(linebuf) >= 4:
            yield (linebuf[0][1:-1], linebuf[1][:-1], linebuf[3][:-1])
            del linebuf[:]

