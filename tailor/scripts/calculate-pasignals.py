#!/usr/bin/env python3
#
# Copyright (c) 2014-5 Institute for Basic Science
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

from tailor.parallel import open_tabix_parallel
from tailor.parsers import parse_sqi
from tailor.sequtils import reverse_complement
from tailor.fileutils import TemporaryDirectory, open_bgzip_writer
from tailor.signalproc import TAILseqSignalProcessor
from tailseqext import PolyALocator
from concurrent import futures
from base64 import b64encode
import pickle
import subprocess as sp
import numpy as np
import pandas as pd
import sys


def run_single_job(input, output, tile):
  try:
    signalconfig, tstretchadj = pickle.load(open(options.paramsfile, 'rb'))
    tstretch_width = (signalconfig['high_probe_scale_range'].stop -
                      signalconfig['high_probe_scale_range'].start)

    sigproc = TAILseqSignalProcessor(signalconfig['cycle_matrix'],
                    signalconfig['color_matrix'], signalconfig['readno'],
                    signalconfig['read_range'].start, signalconfig['read_range'].stop,
                    signalconfig['spot_norm_length'], signalconfig['low_signal_mask'])

    umi_size = signalconfig['spot_norm_length']
    tailscan = options.tailscan

    weights = {'T': options.score_t, 'N': options.score_n, 'A': options.score_acg,
               'C': options.score_acg, 'G': options.score_acg}
    locate_polya = PolyALocator(weights)
    max_term_mod = options.max_term_mod
    min_polya_len = options.min_polya_len
    readcycle_start = signalconfig['read_range'].start

    writer, proc = open_bgzip_writer(output, mode='t')
    processed = skipped_short = skipped_signal_anomaly = 0

    for row in parse_sqi(input()):
        t_stand_out = sigproc.calc_t_stand_out(row.tile, row.seq, row.intensity)
        if t_stand_out is None:
            skipped_signal_anomaly += 1
            continue

        if tailscan: # scan 3' end modifications following poly(A)
            insert_seq = row.seq[readcycle_start + umi_size:]

            seq_start_pa, seq_end_pa = locate_polya(insert_seq, max_term_mod)
            polya_len_seqbased = seq_end_pa - seq_start_pa
            if polya_len_seqbased < min_polya_len:
                # skip calculating PA score for short extremely short poly(A) tails.
                skipped_short += 1
                continue

            adjsignal = t_stand_out[umi_size + seq_start_pa:]
        else: # spike-ins (without terminal modification scanning)
            adjsignal = t_stand_out[umi_size:]

            seq_start_pa, seq_end_pa = locate_polya(insert_seq, max_term_mod)
            seq_start_pa = 0 # fix the starting position of poly(A) for spike-ins
            polya_len_seqbased = seq_end_pa

        if len(adjsignal) >= tstretch_width:
            adjsignal[:tstretch_width] *= tstretchadj[row.tile]
        else:
            adjsignal *= tstretchadj[row.tile][:len(adjsignal)]

        if adjsignal is None or len(adjsignal) < 1:
            skipped_signal_anomaly += 1
            continue

        paencoded = b64encode(adjsignal.astype(np.float32).tostring()).decode('ascii')

        print(row.tile, row.cluster, seq_start_pa, polya_len_seqbased, paencoded,
              sep='\t', file=writer)

        processed += 1

    writer.close()
    proc.wait()

    return (tile, processed, skipped_short, skipped_signal_anomaly)

  except:
    import traceback
    traceback.print_exc()
    raise

def run_jobs(options, output='/dev/stdout'):
    with futures.ProcessPoolExecutor(options.parallel) as executor, \
            TemporaryDirectory(asobj=True) as workdir:
        jobs = []

        tabix_openers = open_tabix_parallel(options.infile[0], named=True)
        for tile, inp in tabix_openers.items():
            joboutput = '{}/{}'.format(workdir.path, tile)
            job = executor.submit(run_single_job, inp, joboutput, tile)
            jobs.append(job)

        countstats = []
        for j in jobs:
            res = j.result()
            countstats.append(res)

        workdir.merge_bgzf(output)

        if options.output_stats is not None:
            statstbl = pd.DataFrame(countstats,
                            columns=['tile', 'processed', 'skipped_short',
                                     'skipped_signal_anomaly']).sort_values(by='tile')
            statstbl.to_csv(options.output_stats, index=False)


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Calculates PA score for intensities')
    parser.add_argument(dest='infile', metavar='FILE', type=str, nargs=1,
                        help='Path to sqi file compressed with bgzip, tabix-indexed')
    parser.add_argument('--parallel', dest='parallel', metavar='N', type=int,
                        help='Number of parallel processes (or threads)', default=4)
    parser.add_argument('--scaling-params', dest='paramsfile', metavar='FILE', type=str,
                        required=True,
                        help='Path to a pickle file that contains per-scale normalization '
                             'matrices in the T-stretch region.')
    parser.add_argument('--no-tail-modification', dest='tailscan', action='store_false',
                        default=True,
                        help='Do not scan poly(A) tail modifications. (for spike-ins)')
    parser.add_argument('--score-T', dest='score_t', metavar='N',
                        default=1, type=int, help='Alignment score for T')
    parser.add_argument('--score-N', dest='score_n', metavar='N',
                        default=-5, type=int, help='Alignment score for N')
    parser.add_argument('--score-ACG', dest='score_acg', metavar='N',
                        default=-10, type=int, help='Alignment score for A, C or G')
    parser.add_argument('--max-terminal-modification', dest='max_term_mod', metavar='N',
                        default=30, type=int,
                        help='Maximum length of 3\' terminal modification after poly(A) tails')
    parser.add_argument('--min-polyA-length', dest='min_polya_len', metavar='N',
                        default=8, type=int,
                        help='Minimum poly(A) length to produce a PA score tuple')
    parser.add_argument('--output-stats', dest='output_stats', metavar='PATH',
                        default=None,
                        help='Write read counting statistics per data flow')

    options = parser.parse_args()

    return options


if __name__ == '__main__':
    options = parse_arguments()
    run_jobs(options)

