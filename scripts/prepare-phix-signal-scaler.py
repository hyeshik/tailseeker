#!/usr/bin/env python3
#
# Copyright (c) 2014-6 Institute for Basic Science
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
from tailseeker.fileutils import TemporaryDirectory
from concurrent import futures
import subprocess as sp
import numpy as np
import pandas as pd
from numpy import linalg
import tables
import os
import pickle
import sys
from itertools import chain
from collections import defaultdict


def load_decrosstalk_matrices(options):
    colormtx = pickle.load(open(options['color_matrix'], 'rb'))

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
            # no positive positions.
            return None
        negmean = adjusted_signals[noncallpositions, basei].mean()

        signal_params.append([negmean, posmean - negmean])

        if posmean - negmean < 1:
            # positive signals are almost indistinguishable from negative signals.
            return None

    return signal_params


def normalize_spot_intensity(intensity, decrosstalk_matrix, norm_params, interval, debug=False):
    adjusted_signals = np.dot(intensity[interval], decrosstalk_matrix)
    return np.array([(adjusted_signals[:, basei] - basalvalue) / dynrange
                     for basei, (basalvalue, dynrange) in enumerate(norm_params)]).T


def open_signal_iterator(inputfilename):
    inputfile = tables.open_file(inputfilename, 'r')
    seqquals = inputfile.get_node('/primary_sqi/PhiX/seqqual')
    intensities = inputfile.get_node('/primary_sqi/PhiX/intensities')

    return inputfile, zip(seqquals.iterrows(), intensities.iterrows())


def collect_signal_stats(options, tileid, inputfile, outputname):
    signalpool = defaultdict(list)
    readnumber = options['readno']
    decrosstalk_matrix = load_decrosstalk_matrices(options)[tileid, readnumber]

    base2i = {b'A'[0]: 0, b'C'[0]: 1, b'G'[0]: 2, b'T'[0]: 3}
    base2iexcluded = {b'A'[0]: [1, 2, 3], b'C'[0]: [0, 2, 3],
                      b'G'[0]: [0, 1, 3], b'T'[0]: [0, 1, 2]}

    read_range = list(options['read_range'])
    read_range[0] -= 1
    read_slice = slice(*read_range)
    spotreference_slice = slice(read_range[0], read_range[0] + options['spot_norm_length'])

    inputtable, itersigs = open_signal_iterator(inputfile)

    for seqqual, intensities in itersigs:
        seq = seqqual['seq']

        norm_params = get_calibration_means(seq, intensities,
                                            decrosstalk_matrix, spotreference_slice)
        if norm_params is None: # no positive signal from any cycles for one or more channels
            continue

        normalized_signal = normalize_spot_intensity(intensities, decrosstalk_matrix,
                                                     norm_params, read_slice)

        for cycleno, (base, normsig) in enumerate(zip(seq[read_slice], normalized_signal)):
            signalpool[base, cycleno].append(normsig.tostring())

    cyclenorm_params = np.nan * np.ones([read_range[1] - read_range[0], 4, 2])
    bad_cycles = []
    num_samples = []

    for cycle in range(read_range[1] - read_range[0]):
        cycleno_full = cycle + read_range[0]
        cyclesignals = []

        for base in b'ACGT':
            signals = np.fromstring(b''.join(signalpool[base, cycle]), np.float64)
            cyclesignals.append(signals.reshape([len(signals) // 4, 4]))

        # Calculate mean signals for positively or negatively called positions for
        # each of A, C, G, and T in each cycle.
        bad_cycle = False
        for basei, base in enumerate(b'ACGT'):
            if len(cyclesignals[basei]) == 0:
                print('No positive signal was found for tile {}, cycle {},'
                      ' base {}'.format(tileid, cycleno_full, base), file=sys.stderr)
                bad_cycle = True
                break
            posmean = cyclesignals[basei][:, basei].mean()
            negmean = np.mean(list(chain(*[cyclesignals[other][:, basei]
                                           for other in base2iexcluded[base]])))

            cyclenorm_params[cycle, basei] = negmean, posmean - negmean

        if bad_cycle:
            bad_cycles.append(cycleno_full)
            num_samples.append(0)
        else:
            num_samples.append(min(map(len, cyclesignals)))

    pickle.dump((tileid, cyclenorm_params, bad_cycles, num_samples), open(outputname, 'wb'))

    inputtable.close()

    return outputname


def write_out_signal_parameters(iter_results, outputfile):
    with tables.open_file(outputfile, 'w') as outf:
        sigprocgrp = outf.create_group('/', 'signalproc',
                                       'Parameters related to signal processing')
        phixgrp = outf.create_group('/signalproc', 'phix',
                                    'Signal scale normalization parameters from PhiX')

        scalegrp = outf.create_group(phixgrp, 'scaling', 'Scaling normalization ranges')
        badcyclesgrp = outf.create_group(phixgrp, 'badcycles',
                                         'List of cycles with low quality signal')
        refcountsgrp = outf.create_group(phixgrp, 'refcounts',
                                         'Number of reference spots')

        # Save this for the CSV output later
        badcycles_ret, refcounts_ret = {}, {}

        for tileid, scaling, badcycles, refcounts in iter_results:
            tilename = 'tile_{}'.format(tileid)

            scaling_arr = outf.create_array(scalegrp, tilename, scaling)
            badcycles_arr = outf.create_array(badcyclesgrp, tilename, badcycles)
            refcounts_arr = outf.create_array(refcountsgrp, tilename, refcounts)

            badcycles_ret[tileid] = badcycles
            refcounts_ret[tileid] = refcounts

    return badcycles_ret, refcounts_ret


def run(options):
    executorclass = (
        futures.ThreadPoolExecutor if options['parallel'] == 1 else futures.ProcessPoolExecutor)

    with executorclass(options['parallel']) as executor, \
            TemporaryDirectory(asobj=True) as tmpdir:
        jobs = []

        for tileid, filename in options['infiles'].items():
            outputname = tmpdir.next_output_file()
            job = executor.submit(collect_signal_stats, options, tileid, filename, outputname)
            jobs.append(job)

        iter_results = (pickle.load(open(j.result(), 'rb')) for j in jobs)
        badcycles, refcounts = write_out_signal_parameters(iter_results, options['output'])

        if options['output_sample_count'] is not None:
            refcountstbl = pd.DataFrame.from_items(sorted(refcounts.items()))
            refcountstbl.index += options['read_range'][0]
            refcountstbl.index.rename('cycle', inplace=True)
            refcountstbl.to_csv(options['output_sample_count'])

        if options['output_bad_cycles'] is not None:
            badcyclestbl = pd.DataFrame.from_items([
                                (tile, [len(cycles), ' '.join(map(str, cycles))])
                                 for tile, cycles in sorted(badcycles.items())]).T
            badcyclestbl.columns = ('count', 'cycle numbers')
            badcyclestbl.index.rename('tile', inplace=True)
            badcyclestbl.to_csv(options['output_bad_cycles'])


if is_snakemake_child:
    options = {
        'infiles': dict(zip(params.tileids, input.signals)),
        'parallel': threads,
        'output': output.paramout,
        'readno': readno,
        'color_matrix': input.colormatrix,
        'read_range': (readstart, readend),
        'spot_norm_length': spotnormlen,
        'output_sample_count': output.refcntstats_out,
        'output_bad_cycles': output.badcycles_out,
    }

    run(options)

