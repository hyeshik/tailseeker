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

from tailseeker.parsers import parse_sam
from tailseeker.tabledefs import *
import subprocess as sp
from Levenshtein import distance as editdist
from collections import deque
import numpy as np
import os

SAMTOOLS_CMD = os.environ.get('SAMTOOLS_CMD', 'samtools')

def parse_sam_options(samline):
    ret = {}
    for token in samline.line[:-1].split(b'\t')[11:]:
        name = token[:2].decode()
        value = int(token[5:].decode()) if token[3] == b'i'[0] else token[5:]
        ret[name] = value
    return ret

def calculate_tag_prority(flags):
    return (
        (0 == (flags & PAFLAG_BARCODE_HAS_MISMATCHES)) * 1024 +
        (0 == (flags & PAFLAG_DELIMITER_NOT_FOUND)) * 512 +
        (0 == (flags & PAFLAG_DELIMITER_HAS_MISMATCH)) * 256 +
        (0 == (flags & PAFLAG_MEASURED_FROM_FLUORESCENCE)) * 128 +
        (0 == (flags & PAFLAG_DARKCYCLE_OVER_THRESHOLD)) * 64 +
        (0 == (flags & PAFLAG_BALANCER_BIASED)) * 32 +
        (0 == (flags & PAFLAG_BALANCER_SIGNAL_BAD)) * 16 +
        (0 == (flags & PAFLAG_DARKCYCLE_EXISTS)) * 8 +
        (0 == (flags & PAFLAG_DELIMITER_IS_SHIFTED)) * 4 +
        (0 == (flags & PAFLAG_NO_POLYA_DETECTED)) * 2 +
        (0 == (flags & PAFLAG_MEASURED_USING_NAIVE_RULER)) * 1)

def reduce_duplicates(group):
    priorities = [calculate_tag_prority(opts['ZF']) for _, opts in group]
    max_priority = max(priorities)

    besttags = [(parsed, opts, prio)
                for (parsed, opts), prio in zip(group, priorities)
                if prio == max_priority]
    besttags.sort(key=lambda vk: vk[1]['Za'])

    polya_lengths = np.array([tag[1]['Za'] for tag in besttags])
    num_clones = np.array([tag[1]['ZD'] for tag in besttags])

    mean_polya_len = (num_clones * polya_lengths).sum() / num_clones.sum()
    nearest = np.abs(polya_lengths - mean_polya_len).argmin()

    total_clones = sum(opts['ZD'] for _, opts in group)

    return (besttags[nearest][0].qname, int(round(mean_polya_len)), total_clones)

def main(options):
    tolerated_interval = options.tolerated_interval
    tolerated_editdist = options.tolerated_editdist
    samproc = sp.Popen([SAMTOOLS_CMD, 'view', '-F', '4', options.bam], stdout=sp.PIPE)

    readstosee = deque()

    def flush():
        reads = readstosee.popleft()
        if len(reads) > 1:
            selected_qname, polya_len, clones = reduce_duplicates(reads)
        else:
            selected_qname = reads[0][0].qname
            polya_len = reads[0][1]['Za']
            clones = reads[0][1]['ZD']

        print(polya_len, clones, selected_qname.decode(), sep='\t')

    for line in parse_sam(samproc.stdout):
        while readstosee:
            first = readstosee[0][0][0]
            if first.rname != line.rname or first.pos + tolerated_interval < line.pos:
                flush()
            else:
                break

        lineopts = parse_sam_options(line)
        for group in readstosee:
            grpfingerprint = group[0][1]['ZM']
            fpdist = editdist(grpfingerprint, lineopts['ZM'])
            if fpdist <= tolerated_editdist:
                group.append((line, lineopts))
                break
        else:
            readstosee.append([(line, lineopts)])

    while readstosee:
        flush()

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description=
                            'Detect duplicated reads by error-tolerated measures')
    parser.add_argument('--bam', dest='bam', type=str,
                        required=True, help='Path to a BAM file')
    parser.add_argument('--tolerated-interval', dest='tolerated_interval',
                        type=int, default=4,
                        help='Maximum coordinate distance inside a clone group')
    parser.add_argument('--tolerated-edit-distance', dest='tolerated_editdist',
                        type=int, default=3,
                        help='Maximum edit distance to the initiator of a clone group')

    return parser.parse_args()


if __name__ == '__main__':
    options = parse_arguments()
    main(options)

