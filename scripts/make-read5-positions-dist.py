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

import subprocess as sp
from tailseeker.parsers import parse_sam
from itertools import groupby
from collections import defaultdict
import pandas as pd

transcript_size = {
    key.encode(): int(value)
    for key, value in (
        line.split() for line
        in open(snakemake.input.transcript_sizes))
}

def iterparse_positions(filename):
    readerproc = sp.Popen(['samtools', 'view', '-F', '20', filename], stdout=sp.PIPE)
    for row in parse_sam(readerproc.stdout):
        trsize = transcript_size[row.rname]
        dist = trsize - row.pos
        yield row.qname, dist

distin = iterparse_positions(snakemake.input.bam)
dcnt = defaultdict(int)
for qname, dists in groupby(distin, key=lambda x: x[0]):
    mindist = min(d for _, d in dists)
    dcnt[mindist] += 1

pd.DataFrame({'count': pd.Series(dcnt)}).to_csv(snakemake.output[0], header=False)

