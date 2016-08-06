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

from tailseeker.fileutils import ParsedLineComment
from tailseeker.parsers import parse_sam
import subprocess as sp
from collections import deque
import numpy as np
import os
import sys

F_UNMAPPED = 4

SAMTOOLS_CMD = os.environ.get('SAMTOOLS_CMD', 'samtools')

def parse_sam_options(samline):
    ret = {}
    for token in samline.line[:-1].split(b'\t')[11:]:
        name = token[:2].decode()
        value = int(token[5:].decode()) if token[3] == b'i'[0] else token[5:]
        ret[name] = value
    return ret

def load_duplicates(duplicates_file):
    r = {}
    for line in open(duplicates_file):
        polya_len_1, polya_len_2, clones, readid = line.split()
        r[readid.encode()] = int(polya_len_1), int(polya_len_2), int(clones)
    return r

def main(options):
    samproc = sp.Popen([SAMTOOLS_CMD, 'view', '-h', options.bam], stdout=sp.PIPE)
    duplicates = load_duplicates(options.duplicates_file)
    output = os.fdopen(sys.stdout.fileno(), 'wb')

    tags_being_processed = [b'ZA', b'Za', b'ZD']

    for row in parse_sam(samproc.stdout):
        if isinstance(row, ParsedLineComment) or row.flag & F_UNMAPPED:
            output.write(row.line)
            # pass unmapped reads and header lines through
            continue
        elif row.qname not in duplicates:
            # skip merged reads
            continue

        dupinfo = duplicates[row.qname]
        if dupinfo[1] == 1:
            # single clone, no modification to tag is required.
            output.write(row.line)
            continue

        samfields = row.line[:-1].split(b'\t')
        samfields = [f for f in samfields if f[:2] not in tags_being_processed]
        samfields.append(b'ZA:i:' + str(dupinfo[0]).encode())
        samfields.append(b'Za:i:' + str(dupinfo[1]).encode())
        samfields.append(b'ZD:i:' + str(dupinfo[2]).encode() + b'\n')
        output.write(b'\t'.join(samfields))

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description=
                            'Filter detected duplicated reads by error-tolerated measures')
    parser.add_argument('--bam', dest='bam', type=str,
                        required=True, help='Path to a BAM file')
    parser.add_argument('--duplicates', dest='duplicates_file',
                        required=True, help='Path to a table of detected duplicates')

    return parser.parse_args()


if __name__ == '__main__':
    options = parse_arguments()
    main(options)

