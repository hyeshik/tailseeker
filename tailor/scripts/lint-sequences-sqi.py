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

from tailor.fileutils import (
    TemporaryDirectory, open_gzip_buffered, MultiJoinIterator, open_bgzip_writer)
from tailor.parallel import open_tabix_parallel, TabixOpener
from tailor.parsers import parse_sqi_lite
from collections import Counter
from concurrent import futures
from operator import itemgetter
import subprocess as sp
import os


def string_dist(s1, s2):
    return sum(a != b for a, b in zip(s1, s2))


def parse_sqi_fast(opener):
    for line in opener():
        tile, cluster, _ = line.split(b'\t', 2)
        yield (tile, int(cluster)), line

def parse_idlist(opener):
    for line in opener():
        tile, cluster = line.split(b'\t')
        yield (tile, int(cluster))

def make_it_lint(options, input, outputname):
  try:
    output = open_bgzip_writer(outputname, 'b')

    sqi_in = parse_sqi_fast(input)
    sqi_key = itemgetter(0)

    whitelist_in = parse_idlist(TabixOpener(options.idlist, input.interval))
    whitelist_key = lambda x: x

    # bind preamble filtering-related options into pre-calculated local variables
    if options.preamblepos is not None:
        preamblepos, preambleseq = options.preamblepos - 1, options.preambleseq
        preambleslice = slice(preamblepos, preamblepos + len(preambleseq))
    else:
        preamblepos = preambleseq = None
    preamblemismatch = options.preamblemismatch
    preambleend= options.preambleend # use only if all the other related variables are not defined

    # bind balance region check options into pre-calculated local variables
    if options.balanceregion is not None:
        balance_start1, balance_end = map(int, options.balanceregion.split(':'))
        balanceslice = slice(balance_start1 - 1, balance_end)
        balancemin = options.balancemin
    else:
        balanceslice = balancemin = None

    for key, sqiit, whitelistit in MultiJoinIterator([sqi_in, whitelist_in],
                                                     [sqi_key, whitelist_key]):
        sqiline = next(sqiit, None)
        whitelistline = next(whitelistit, None)
        if sqiline is None or whitelistline is None:
            continue

        tile, cluster = key
        line = sqiline[1]

        spot = next(parse_sqi_lite([line]))

        # check for preamble sequences
        if preamblepos is not None:
            if string_dist(spot.seq[preambleslice], preambleseq) > preamblemismatch:
                continue

        # check signal balance control region
        if balanceslice is not None:
            bseq = spot.seq[balanceslice]
            if min(bseq.count(c) for c in 'ACGT') < balancemin:
                continue

        output.write(line)

    output.close()

    return outputname
  except:
    import traceback
    traceback.print_exc()
    raise


def run(options):
    infile = options.infile[0]

    with futures.ProcessPoolExecutor(options.parallel) as executor, \
            TemporaryDirectory(asobj=True) as tmpdir:
        jobs = []

        for jobno, opener in enumerate(open_tabix_parallel(infile)):
            outputname = tmpdir.next_output_file()
            job = executor.submit(make_it_lint, options, opener, outputname)
            jobs.append(job)

        for j in jobs:
            j.result()

        tmpdir.merge_bgzf(options.output)


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Filters duplicated reads from sqi')
    parser.add_argument(dest='infile', metavar='FILE', type=str, nargs=1,
                        help='Path to a sqi file')
    parser.add_argument('--id-list', dest='idlist', metavar='FILE', type=str,
                        required=True, help='Cluster names to keep')
    parser.add_argument('--preamble-position', dest='preamblepos', metavar='N', type=int,
                        default=None, help='1-based position of preamble sequence')
    parser.add_argument('--preamble-sequence', dest='preambleseq', metavar='SEQ', type=str,
                        default=None, help='Sequence string to check for preamble')
    parser.add_argument('--preamble-mismatch', dest='preamblemismatch', metavar='N', type=int,
                        default=1, help='Maximum mismatches to filter by preamble sequence')
    parser.add_argument('--balance-region', dest='balanceregion', metavar='start:end',
                        type=str, default=None,
                        help='1-based coordinate for interval for balance check')
    parser.add_argument('--balance-minimum', dest='balancemin', metavar='N', type=int,
                        default=None, help='Minimum count to filter in balance region')
    parser.add_argument('--preamble-end', dest='preambleend', metavar='N', type=int,
                        default=None, help='End position (1-based) of preamble')
    parser.add_argument('--parallel', dest='parallel', metavar='N', type=int,
                        help='Number of parallel processors', default=4)
    parser.add_argument('--output', dest='output', metavar='FILE', type=str,
                        help='Path to output file (bgz)', required=True)
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    options = parse_arguments()
    run(options)

