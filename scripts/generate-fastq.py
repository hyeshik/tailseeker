#!/usr/bin/env python3
#
# Copyright (c) 2015 Institute for Basic Science
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

from tailor.parsers import parse_sqi_lite, parse_polya_calls
from tailor.parallel import open_tabix_parallel
from tailor.fileutils import TemporaryDirectory, MultiJoinIterator
from tailseqext import SimpleBGZFWriter
from concurrent import futures
from collections import Counter
import sys


def run_subjob(sqiopener, paopener, reads, joboutput):
  try:
    sqi_in = parse_sqi_lite(sqiopener())
    pa_in = parse_polya_calls(paopener())
    common_key = lambda s: (s.tile, s.cluster)

    preader = MultiJoinIterator([sqi_in, pa_in], [common_key, common_key])
    reads = [
        [readspec['name'], slice(readspec['start'], readspec['end']),
         readspec['trim'],
         SimpleBGZFWriter(joboutput.format(readname=readspec['name']))]
        for readspec in reads
    ]

    for (tile, cluster), sqispot, paspot in preader:
        sqispot, paspot = next(sqispot), next(paspot)
        coded_name = '{}:{:08d}:{:03d}:{:03d}\n'.format(tile, cluster, paspot.polya_len,
                                                        paspot.start_pos).encode('ascii')

        for name, span, trim, writer in reads:
            if trim is not None:
                delimiterend = sqispot.istart
                span = slice(delimiterend, min(span.stop, delimiterend + trim))

            seq = sqispot.seq[span].encode('ascii')
            qual = sqispot.qual[span]
            writer(b'@' + coded_name)
            writer(seq)
            writer(b'\n+' + coded_name)
            writer(qual)
            writer(b'\n')

    for _, _, _, writer in reads:
        writer.close()

  except:
    import traceback
    traceback.print_exc()
    raise


def terminated_iterator():
    return
    yield None


def run(options, reads):
    sqi_openers = open_tabix_parallel(options.input_sqi, named=True)
    pacall_openers = open_tabix_parallel(options.input_pa_call, named=True)

    tiles = sorted(sqi_openers)
    sqi_openers = [sqi_openers[k] for k in tiles]
    pacall_openers = [(pacall_openers[k] if k in pacall_openers
                       else terminated_iterator) for k in tiles]

    with futures.ProcessPoolExecutor(options.parallel) as executor, \
            TemporaryDirectory(asobj=True) as workdir:
        jobs = []

        for inputno, (sqiopen, paopen) in enumerate(zip(sqi_openers, pacall_openers)):
            joboutput = '{}/{{readname}}_{:08x}'.format(workdir.path, inputno)
            job = executor.submit(run_subjob, sqiopen, paopen, reads, joboutput)
            jobs.append(job)

        for j in jobs:
            j.result()

        for readspec in reads:
            workdir.merge_bgzf(options.output.replace('XX', readspec['name']),
                               prefix='{}_'.format(readspec['name']))


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Produce FASTQ sequence files with measured '
                                                 'poly(A) tail information')
    parser.add_argument('readspan', type=str, nargs='+',
                        help='Ranges of reads, specify in name,start,end,trim (1-based).'
                        ' (trim: N=don\'t use delimiter, number=number of cycles to keep)')
    parser.add_argument('--input-sqi', dest='input_sqi', metavar='FILE', required=True,
                        type=str, help='Path to a file containing sequence and intensity')
    parser.add_argument('--input-pa-call', dest='input_pa_call', metavar='FILE', required=True,
                        type=str, help='Path to a file containing poly(A) calls')
    parser.add_argument('--parallel', dest='parallel', metavar='N',
                        default=8, type=int, help='Number of parallel processes')
    parser.add_argument('--output', dest='output', metavar='FILE', type=str,
                        required=True, help='Path to output file (fastq.gz format), '
                                            'include XX for read number.')

    options = parser.parse_args()

    reads = []
    for readspan in options.readspan:
        fields = readspan.split(',')
        if len(fields) != 4:
            print('Reads must be specified in name,start,end,[yn] format.', file=sys.stderr)
            sys.exit(1)

        reads.append({'name': fields[0],
                      'start': int(fields[1]) - 1,
                      'end': int(fields[2]),
                      'trim': None if fields[3] == 'N' else int(fields[3])})

    return options, reads


if __name__ == '__main__':
    options, reads = parse_arguments()
    run(options, reads)

