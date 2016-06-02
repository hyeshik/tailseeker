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

from tailseeker.parsers import parse_sam, parse_refined_taginfo
from tailseeker.fileutils import MultiJoinIterator, ParsedLineComment
import subprocess as sp
import sys
import gzip
import os

"""
 48 parse_refined_taginfo = LineParser([
 49     ('tile', None),
 50     ('cluster', int),
 51     ('pflags', int),
 52     ('clones', int),
 53     ('polyA', int),
 54     ('U', int),
 55     ('G', int),
 56     ('C', int),
 57     ('mods', None),
 58 ], linefeed=b'\n')
"""


def process(taginfo_file):
    taginfo_input = gzip.open(taginfo_file, 'rb')
    sam_input = os.fdopen(sys.stdin.fileno(), 'rb')

    # parsing iterators
    samit = parse_sam(sam_input)
    taginfoit = parse_refined_taginfo(taginfo_input)

    # key functions for joiner
    def parse_readid_from_sam(samrow):
        if isinstance(samrow, ParsedLineComment):
            return b'', 0
        qtokens = samrow.qname.split(b':', 2)
        return (qtokens[0], int(qtokens[1]))
    taginfokey = lambda x: (x.tile, x.cluster)

    joined_it = MultiJoinIterator([samit, taginfoit], [parse_readid_from_sam, taginfokey])
    output = os.fdopen(sys.stdout.fileno(), 'wb')
    tagformat = '\tZF:i:{}\tZD:i:{}\tZU:i:{}\tZG:i:{}\tZC:i:{}\tZa:i:{}\n'

    for (tile, cluster), samrows, taginforows in joined_it:
        samrows = list(samrows)
        taginforows = list(taginforows)

        if len(taginforows) < 1:
            output.write(b''.join(row.line for row in samrows))
            continue
        elif len(samrows) < 1:
            continue

        taginfo = taginforows[0]
        tags_formatted = tagformat.format(taginfo.pflags, taginfo.clones, taginfo.U,
                                          taginfo.G, taginfo.C, taginfo.polyA)

        for row in samrows:
            output.write(row.line[:-1])
            output.write(tags_formatted.encode())


if __name__ == '__main__':
    taginfo = sys.argv[1]
    process(taginfo)

