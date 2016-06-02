#!/usr/bin/env python3
#
# Copyright (c) 2013-2015 Institute for Basic Science
# Copyright (c) 2011-2013 Hyeshik Chang
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

import re
import os.path

__all__ = [
    'reverse_complement',
    'reverse_complement_bytes',
    'GiantFASTAFile',
]

revcmptrans = str.maketrans('ATUGCatugc', 'TAACGtaacg')
def reverse_complement(seq):
    return seq.translate(revcmptrans)[::-1]

revcmptrans_bytes = bytes.maketrans(b'ATUGCatugc', b'TAACGtaacg')
def reverse_complement_bytes(seq):
    return seq.translate(revcmptrans_bytes)[::-1]


whitespace = re.compile('[ \t\r\n]')
class GiantFASTAFile(object):

    def __init__(self, filename):
        if not os.path.exists(filename + '.fai'):
            import pysam
            pysam.faidx(filename)

        self.fasta = open(filename)
        self.index = self.load_index(filename + '.fai')

    def load_index(self, filename):
        index = {}
        for line in open(filename):
            fields = line[:-1].split('\t')
            index[fields[0]] = tuple(map(int, fields[1:]))
        return index

    def get(self, seqid, start=None, stop=None, strand='+'): # zero-based, half-open
        length, filepos, colwidth, linesize = self.index[seqid]

        if start is None and stop is None:
            offset_st = filepos
            linenum_en = length // colwidth
            offset_en = filepos + length + linenum_en * (linesize - colwidth)
        else:
            start = max(0, start)
            stop = min(length, stop)

            linenum_st = start // colwidth
            offset_st = filepos + start + linenum_st * (linesize - colwidth)
            linenum_en = stop // colwidth
            offset_en = filepos + stop + linenum_en * (linesize - colwidth)

        self.fasta.seek(offset_st, 0)
        seq = whitespace.sub('', self.fasta.read(offset_en - offset_st))
        return seq if strand == '+' else reverse_complement(seq)
