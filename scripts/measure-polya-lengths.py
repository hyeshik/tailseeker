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

from tailseeker.powersnake import *
from snakemake.shell import shell
from snakemake.utils import format
from tailseeker.fileutils import TemporaryDirectory
import numpy as np
import struct
import gzip

BINDIR = params.BINDIR
BGZIP_CMD = params.BGZIP_CMD
CONF = params.CONF

def load_sig_dists(filename):
    inpf = gzip.open(filename, 'rb')
    elemsize, cycles, bins = struct.unpack('<III', inpf.read(12))
    assert elemsize == 4

    counts = np.frombuffer(inpf.read(), np.uint32).reshape((cycles, bins))
    return counts.astype(np.uint64)

def write_sig_dists(counts, filename):
    with gzip.open(filename, 'wb') as outf:
        outf.write(struct.pack('<III', 4, counts.shape[0], counts.shape[1]))
        outf.write(counts.astype(np.uint32).tostring())

counts_aggr = None

for signals, taginfo, out in zip(input.signals, input.taginfo,
                                 output.taginfo):
    cmd = format('{BINDIR}/tailseq-polya-ruler {wildcards.tile} {signals} \
        {input.score_cutoffs} {CONF[polyA_finder][signal_analysis_trigger]} \
        {CONF[polyA_ruler][downhill_extension_weight]} \
        {taginfo} {CONF[polyA_seeder][dist_sampling_bins]} \
        {CONF[polyA_ruler][signal_resampling_gap]} \
        {output.sigdists} | {BGZIP_CMD} -c > {out}', wildcards=wildcards,
        input=input, output=output)
    shell(cmd)

    counts_new = load_sig_dists(output.sigdists)
    if counts_aggr is None:
        counts_aggr = counts_new
    else:
        counts_aggr += counts_new

write_sig_dists(counts_aggr, output.sigdists)

