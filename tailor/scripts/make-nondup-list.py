#!/usr/bin/env python3
#
# Copyright (c) 2013-5 Institute for Basic Science
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

from tailor.fileutils import open_gzip_buffered, MultiJoinIterator
import subprocess as sp
from io import TextIOWrapper
import os

def unpack_readid(s):
    tile, cluster = s.split()
    return (tile, int(cluster))

def open_blacklist(filename):
    return (unpack_readid(line) for line in open_gzip_buffered(filename, 'rt'))

def open_fastq_idonly(filename):
    cmd = """zcat {input} | awk -F: '{{ if (NR % 4 == 1) print substr($1, 2), 0+$2; }}'"""
    proc = sp.Popen(cmd.format(input=filename), shell=True, stdout=sp.PIPE)
    for line in TextIOWrapper(proc.stdout):
        yield unpack_readid(line)

def run(options):
    fastqids_in = open_fastq_idonly(options.fastq)
    excludes_in = [open_blacklist(excfile) for excfile in options.exclude]

    for grps in MultiJoinIterator([fastqids_in] + excludes_in,
                                  [lambda x: x] * (len(excludes_in) + 1)):
        fqkey = next(grps[1], None)
        if fqkey is None or any(next(g, None) is not None for g in grps[2:]):
            continue

        print(*fqkey, sep='\t')


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Make non-duplicated list of reads')
    parser.add_argument('--fastq', dest='fastq', metavar='FILE', type=str,
                        required=True, help='Path to a fastq file with all read 1 sequences')
    parser.add_argument('--exclude', dest='exclude', metavar='FILE', action='append',
                        default=[], help='Path to blacklist of seq IDs')
    options = parser.parse_args()

    return options

if __name__ == '__main__':
    options = parse_arguments()
    run(options)
