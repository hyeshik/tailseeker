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

from tailor.fileutils import TemporaryDirectory
from tailor.parsers import parse_sqi
from tailor.parallel import open_tabix_parallel, TabixOpener
from concurrent import futures
import subprocess as sp
import numpy as np
from numpy import linalg
import os
import glob
import pickle
import sys
import csv
import re
from itertools import chain
from collections import defaultdict


def load_decrosstalk_matrices(options):
    colormtx = pickle.load(open(options.colormtx, 'rb'))

    decrosstalk_matrices = {}
    for (vtile, readno), mtx in colormtx.items():
        decrosstalk_matrices[vtile, int(readno)] = linalg.inv(mtx.T)

    return decrosstalk_matrices


def get_calibration_means(seq, intensity, decrosstalkmtx, interval):
    seqarr = np.fromstring(seq[interval], dtype=np.uint8)
    adjusted_signals = np.dot(intensity[interval], decrosstalkmtx)

    signal_params = []

    for basei, base in enumerate(b'ACGT'):
        callpositions = np.where(seqarr == base)[0]
        if len(callpositions) < 1:
            return None

        posmean = adjusted_signals[callpositions, basei].mean()

        noncallpositions = np.where(seqarr != base)[0]
        if len(noncallpositions) < 1:
            return None
        negmean = adjusted_signals[noncallpositions, basei].mean()

        signal_params.append([negmean, posmean - negmean])

    return signal_params


def normalize_spot_intensity(intensity, decrosstalk_matrix, norm_params, interval):
    adjusted_signals = np.dot(intensity[interval], decrosstalk_matrix)
    return np.array([(adjusted_signals[:, basei] - basalvalue) / dynrange
                     for basei, (basalvalue, dynrange) in enumerate(norm_params)]).T


def collect_signal_stats(options, opener, outputname):
    tilenumber = None
    signalpool = defaultdict(list)
    readnumber = options.readno
    decrosstalk_matrices = load_decrosstalk_matrices(options)

    base2i = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    base2iexcluded = {'A': [1, 2, 3], 'C': [0, 2, 3], 'G': [0, 1, 3], 'T': [0, 1, 2]}

    read_range = list(map(int, options.readrange.split(':')))
    read_range[0] -= 1
    read_slice = slice(*read_range)
    spotreference_slice = slice(read_range[0], read_range[0] + options.spotnormlen)

    for spot in parse_sqi(opener()):
        if tilenumber is None:
            tilenumber = spot.tile
            decrosstalk_matrix = decrosstalk_matrices[tilenumber, readnumber]

        norm_params = get_calibration_means(spot.seq, spot.intensity,
                                            decrosstalk_matrix, spotreference_slice)
        if norm_params is None: # no positive signal from any cycles for one or more channels
            continue

        normalized_signal = normalize_spot_intensity(spot.intensity, decrosstalk_matrix,
                                                     norm_params, read_slice)

        for cycleno, (base, normsig) in enumerate(zip(spot.seq[read_slice],
                                                      normalized_signal)):
            signalpool[base, cycleno].append(normsig.tostring())


    cyclenorm_params = {}
    bad_cycles = []
    num_samples = []

    for cycle in range(read_range[1] - read_range[0]):
        cycleno_full = cycle + read_range[0]
        cyclesignals = []

        for base in 'ACGT':
            signals = np.fromstring(b''.join(signalpool[base, cycle]), np.float64)
            cyclesignals.append(signals.reshape([len(signals) // 4, 4]))

        # Calculate mean signals for positively or negatively called positions for
        # each of A, C, G, and T in each cycle.
        bad_cycle = False
        for basei, base in enumerate('ACGT'):
            if len(cyclesignals[basei]) == 0:
                print('No positive signal was found for tile {}, cycle {},'
                      ' base {}'.format(tilenumber, cycle, base), file=sys.stderr)
                bad_cycle = True
                break
            posmean = cyclesignals[basei][:, basei].mean()
            negmean = np.mean(list(chain(*[cyclesignals[other][:, basei]
                                           for other in base2iexcluded[base]])))

            cyclenorm_params[cycleno_full, base] = negmean, posmean - negmean

        if bad_cycle:
            bad_cycles.append(cycleno_full)
            num_samples.append(0)
        else:
            num_samples.append(min(map(len, cyclesignals)))


    pickle.dump(({tilenumber: cyclenorm_params},
                 {tilenumber: bad_cycles},
                 {tilenumber: num_samples}), open(outputname, 'wb'))

    return outputname

def run(options):
    infile = options.infile[0]
    #collect_signal_stats(options, open_tabix_parallel(infile)[0], options.output)
    #raise SystemExit

    with futures.ProcessPoolExecutor(options.parallel) as executor, \
            TemporaryDirectory(asobj=True) as tmpdir:
        jobs = []

        for jobno, opener in enumerate(open_tabix_parallel(infile)):
            outputname = tmpdir.next_output_file()
            job = executor.submit(collect_signal_stats, options, opener, outputname)
            jobs.append(job)

        scale_parameter_outputs = {}
        bad_cycles_outputs = {}
        num_samples_outputs = {}
        for j in jobs:
            scaleparams, badcycles, numsamples = pickle.load(open(j.result(), 'rb'))
            scale_parameter_outputs.update(scaleparams)
            bad_cycles_outputs.update(badcycles)
            num_samples_outputs.update(numsamples)

        pickle.dump(scale_parameter_outputs, open(options.output, 'wb'))

        if options.samplecountcsv is not None:
            sc_out = csv.writer(open(options.samplecountcsv, 'w'))
            read_range = list(map(int, options.readrange.split(':')))

            sc_out.writerow(['Tile'] + list(range(read_range[0], read_range[1] + 1)))
            for tile, samplecounts in sorted(num_samples_outputs.items()):
                sc_out.writerow([tile] + samplecounts)


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Calculates normalization parameters '
                                                 'for fluorescence signals')
    parser.add_argument(dest='infile', metavar='FILE', type=str, nargs=1,
                        help='Path to a sqi file')
    parser.add_argument('--parallel', dest='parallel', metavar='NUMBER', type=int, default=1,
                        help='Number of parallel processes')
    parser.add_argument('--output', dest='output', metavar='FILE', type=str, required=True,
                        help='Output file for collected statistics')
    parser.add_argument('--read', dest='readno', metavar='NUMBER', type=int, required=True,
                        help='Read number to process (3 in the standard paired-end setup).')
    parser.add_argument('--color-matrix', dest='colormtx', metavar='FILE', type=str,
                        required=True,
                        help='Path to a pickle file that contains color matrices.')
    parser.add_argument('--read-range', dest='readrange', metavar='start:end', type=str,
                        required=True,
                        help='Cycle number range of the read in 1-based inclusive coordinate.')
    parser.add_argument('--spot-norm-length', dest='spotnormlen', metavar='LEN', type=int,
                        default=20,
                        help='Use the first N cycles for detection of spot brightness')
    parser.add_argument('--sample-number-stats', dest='samplecountcsv', metavar='FILE',
                        type=str, default=None,
                        help='Statistical output for sample counts used in distribution '
                             'estimations.')
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    options = parse_arguments()
    run(options)

