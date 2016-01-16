#!/usr/bin/env python3
#
# Copyright (c) 2013-2015 Institute for Basic Science
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

from tailor.fileutils import TemporaryDirectory
from tailor.parsers import parse_sqi_lite
from tailor.parallel import open_tabix_parallel
import os
import sys
import csv
import subprocess as sp
from itertools import groupby
from concurrent import futures
from functools import partial
from collections import defaultdict
from random import randrange, seed
from io import TextIOWrapper

BGZIP_CMD = 'bgzip'

def generate_id_sequences(input, num, slices, tmpdir):
    outputfile = os.path.join(tmpdir, str(num))
    output = open(outputfile, 'w')

    for spot in parse_sqi_lite(input()):
        seq = spot.seq
        idseq = ''.join(seq[s] for s in slices)
        print(idseq, spot.tile, spot.cluster, sep='\t', file=output)

    output.close()

    return outputfile

def pack_readid(tile, cluster):
    return '{}:{:08d}'.format(tile, int(cluster))


def sort_and_find_duplicates(output, idfiles, parallel, debug_output=None,
                             dupcount_output=None, clonelist_output=None):
    soproc = sp.Popen(['sort', '--parallel={}'.format(parallel), '-k1,1'] + idfiles,
                      stdout=sp.PIPE)
    grpstats = defaultdict(int)

    seed(8809)

    for seq, grp in groupby((line[:-1].split() for line in TextIOWrapper(soproc.stdout)),
                            key=lambda x: x[0]):
        grp = [(tile, cluster) for _, tile, cluster in grp]
        grpstats[len(grp)] += 1
        duplicity = len(grp)

        if duplicity >= 2:
            survivor = randrange(duplicity)
            surviving_cluster = grp.pop(survivor)
            print('\n'.join('{}\t{}'.format(*p) for p in grp), file=output)

            if debug_output is not None and len(grp) >= 10:
                print('{}\t{}\t{}\t{}'.format(*(duplicity, seq) + surviving_cluster),
                      file=debug_output)
        else:
            surviving_cluster = grp.pop()

        if dupcount_output is not None:
            print('{}\t{}\t{}'.format(*(surviving_cluster + (duplicity,))),
                  file=dupcount_output)

        if clonelist_output is not None:
            print(pack_readid(*surviving_cluster),
                  ' '.join(pack_readid(*args) for args in grp),
                  file=clonelist_output)

    return grpstats


def write_dupstats(output, dupstats):
    w = csv.writer(output)
    w.writerow(('duplicated reads', 'groups of duplicated reads'))
    w.writerows(sorted(dupstats.items(), key=lambda x: x[0]))


def run(options, regions):
    infile = options.infile[0]
    id_slices = [slice(begin-1, end) for begin, end in regions]

    if options.dupcounts is not None:
        dupcount_out = os.popen('sort -k1,1 -k2,2n | {bgzip} -c /dev/stdin > {output}'.format(
                                    output=options.dupcounts, bgzip=BGZIP_CMD), 'w')
    else:
        dupcount_out = None

    clonelist_out = open(options.clonelist, 'w') if options.clonelist is not None else None

    with TemporaryDirectory() as tmpdir:
        with futures.ProcessPoolExecutor(options.parallel) as executor:
            jobs = []

            for jobid, opener in enumerate(open_tabix_parallel(infile)):
                job = executor.submit(generate_id_sequences, opener, jobid, id_slices, tmpdir)
                jobs.append(job)

            result_files = [job.result() for job in jobs]

        debug_output = open(options.debug, 'w') if options.debug is not None else None

        dupstats = sort_and_find_duplicates(sys.stdout, result_files, options.parallel,
                                            debug_output, dupcount_out, clonelist_out)

        if options.dupstats is not None:
            write_dupstats(open(options.dupstats, 'w'), dupstats)

    if dupcount_out is not None:
        dupcount_out.close()


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Prepares a death list for duplicated reads')
    parser.add_argument(dest='infile', metavar='FILE', type=str, nargs=1,
                        help='Path to a sqi file')
    parser.add_argument('--region', dest='regions', metavar='start:end', type=str,
                        default=[], action='append', help='Cycle range to be regarded')
    parser.add_argument('--output-stats', dest='dupstats', metavar='FILE', type=str,
                        help='Path to file for duplication statistics (CSV)', default=None)
    parser.add_argument('--parallel', dest='parallel', metavar='N', type=int,
                        help='Number of parallel processors', default=4)
    parser.add_argument('--output-dupcounts', dest='dupcounts', metavar='FILE', type=str,
                        default=None, help='Path to file for counts of duplicated reads')
    parser.add_argument('--output-clone-list', dest='clonelist', metavar='FILE', type=str,
                        default=None, help='Path to file for clones list of duplicated reads')
    parser.add_argument('--output-debug', dest='debug', metavar='FILE', type=str,
                        default=None, help='Path to file for debug output')
    options = parser.parse_args()

    regions = [map(int, reg.split(':')) for reg in options.regions]

    return options, regions


if __name__ == '__main__':
    if 'BGZIP_CMD' in os.environ:
        BGZIP_CMD = os.environ['BGZIP_CMD']

    options, regions = parse_arguments()
    run(options, regions)

