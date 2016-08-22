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
from tailseeker.stats import gaussian_kde
from scipy import optimize
from concurrent import futures
import sys
import numpy as np
import struct
import gzip

NUM_CLOSEST_CYCLES_INTERPOLATION    = 3
NUM_CLOSEST_CYCLES_EXTRAPOLATION    = 5

def load_sig_dists(filename):
    inpf = gzip.open(filename, 'rb')
    elemsize, cycles, bins = struct.unpack('<III', inpf.read(12))
    assert elemsize == 4

    counts = np.fromstring(inpf.read(), np.uint32).reshape((cycles, bins))
    spots_per_cycle = counts.sum(axis=1)
    return counts, spots_per_cycle

def load_all_signal_dists():
    sigdists = {}
    spotcounts = {}
    for f in input:
        tilename = f.rsplit('_', 1)[1].split('.', 1)[0]
        signaltype = f.rsplit('/', 1)[1].split('_', 1)[0]
        sigdists[tilename, signaltype], spotcounts[tilename, signaltype] = \
            load_sig_dists(f)
    return sigdists

def find_eqodds_point(poscounts, negcounts, xsamples):
    poskde = gaussian_kde(xsamples, weights=poscounts, bw_method=KDE_BANDWIDTH)
    negkde = gaussian_kde(xsamples, weights=negcounts, bw_method=KDE_BANDWIDTH)

    logLRfun = lambda x: np.log(max(VERY_SMALL_PROBABILITY, poskde(x)) /
                                max(VERY_SMALL_PROBABILITY, negkde(x)))

    try:
        return optimize.bisect(logLRfun, CUTOFF_RANGE_LOW, CUTOFF_RANGE_HIGH)
    except ValueError:
        return

def find_all_eqodds_point(poscounts, negcounts, xsamples):
    cycles_to_try = np.where((poscounts.sum(axis=1) >= MIN_SPOTS_FOR_DIST) &
                             (negcounts.sum(axis=1) >= MIN_SPOTS_FOR_DIST))[0].tolist()
    return {
        cycle: find_eqodds_point(poscounts[cycle], negcounts[cycle], xsamples)
        for cycle in cycles_to_try}

def infer_missing_cutoffs(cutoffs, allcycles):
    available_cycles = np.array(list(sorted(cutoffs)))
    missing_cycles = np.array(list(sorted(set(allcycles) - set(cutoffs))))
    group_breaks = [0] + (np.where(np.diff(missing_cycles) > 1)[0] + 1).tolist() + [
                    len(missing_cycles)]

    inferred_cutoffs = {}
    for smallest, largest in zip(group_breaks[:-1], group_breaks[1:]):
        cycle_left_inner = missing_cycles[smallest]
        cycle_right_inner = missing_cycles[largest - 1]
        icycles = list(range(cycle_left_inner, cycle_right_inner + 1))

        cycles_on_left = available_cycles[available_cycles < cycle_left_inner].tolist()
        cycles_on_right = available_cycles[available_cycles > cycle_right_inner].tolist()

        if cycles_on_left: # left ok
            if cycles_on_right: # right ok
                # Infer a cutoff with linear interpolation
                leftcycles = cycles_on_left[-NUM_CLOSEST_CYCLES_INTERPOLATION:]
                rightcycles = cycles_on_right[:NUM_CLOSEST_CYCLES_INTERPOLATION]
                leftcutoffs = [cutoffs[c] for c in leftcycles]
                rightcutoffs = [cutoffs[c] for c in rightcycles]
                interpolated_cutoffs = np.interp(
                    icycles, leftcycles + rightcycles, leftcutoffs + rightcutoffs)
                for c, cutoff in zip(icycles, interpolated_cutoffs):
                    inferred_cutoffs[c] = cutoff
            else:
                # Extrapolate from the left-side
                leftcycles = cycles_on_left[-NUM_CLOSEST_CYCLES_EXTRAPOLATION:]
                leftcutoffs = [cutoffs[c] for c in leftcycles]
                extcutoff = np.median(leftcutoffs)
                for c in icycles:
                    inferred_cutoffs[c] = extcutoff
        else:
            # Extrapolate from the right-side
            rightcycles = cycles_on_right[:NUM_CLOSEST_CYCLES_EXTRAPOLATION]
            rightcutoffs = [cutoffs[c] for c in rightcycles]
            extcutoff = np.median(rightcutoffs)
            for c in icycles:
                inferred_cutoffs[c] = extcutoff

    c = cutoffs.copy()
    c.update(inferred_cutoffs)
    return c

def get_report_marks(direct, inferred):
    r = []
    for i in range(TOTAL_CYCLES):
        if i in direct:
            r.append('d')
        elif i in inferred:
            r.append('i')
        else:
            r.append('x')
    return ''.join(r)

def calculate_optimal_cutoffs(sigdists, tiles):
    # Sum the counts from all tiles
    postotal = sum([sigdists[tile['id'], 'pos'] for tile in tiles])
    negtotal = sum([sigdists[tile['id'], 'neg'] for tile in tiles])

    xsamples = np.arange(0, 1, 1 / postotal.shape[1])

    runid = tiles[0]['source']
    all_cycles = np.where(negtotal.sum(axis=1) > 0)[0].tolist()
    runwide_cutoffs = find_all_eqodds_point(postotal, negtotal, xsamples)

    tile_cutoffs = {}
    with futures.ProcessPoolExecutor(NUM_THREADS) as executor:
        for tile in sorted(tiles, key=lambda x: x['id']):
            poscounts = sigdists[tile['id'], 'pos']
            negcounts = sigdists[tile['id'], 'neg']
            job = executor.submit(find_all_eqodds_point, poscounts, negcounts, xsamples)
            job.add_done_callback(lambda f, tileid=tile['id']:
                                    tile_cutoffs.__setitem__(tileid, f.result()))

    runwide_cutoffs_inferred = infer_missing_cutoffs(runwide_cutoffs, all_cycles)
    cutoff_reports = {runid: get_report_marks(runwide_cutoffs, runwide_cutoffs_inferred)}

    tile_cutoffs_inferred = {}
    for tileid, tile_cutoffs in tile_cutoffs.items():
        cutoffs = tile_cutoffs.copy()
        for c in set(all_cycles) - set(tile_cutoffs.keys()):
            cutoffs[c] = runwide_cutoffs_inferred[c]
        tile_cutoffs_inferred[tileid] = cutoffs

        cutoff_reports[tileid] = get_report_marks(tile_cutoffs, cutoffs)

    return tile_cutoffs_inferred, cutoff_reports

def write_outputs(cutoffs, bases):
    with open(output.cutoff_values, 'wt') as outf:
        for tile, cutoffs in sorted(cutoffs.items()):
            print(tile + '\t', end='', file=outf)
            for i in range(TOTAL_CYCLES):
                if i in cutoffs:
                    print(format(cutoffs[i], '.6f') + '\t', end='', file=outf)
                else:
                    print('nan\t', end='', file=outf)
            print(file=outf)

    with open(output.cutoff_bases, 'wt') as outf:
        for tile, basismarks in sorted(bases.items()):
            print(tile, basismarks, sep='\t', file=outf)

def main():
    sigdists = load_all_signal_dists()

    tilelist = {}
    for tileinfo in params.tileinfo.values():
        tilelist.setdefault(tileinfo['source'], [])
        tilelist[tileinfo['source']].append(tileinfo)

    cutoffs, bases = {}, {}

    for source, tiles in tilelist.items():
        tile_cutoffs, tile_bases = calculate_optimal_cutoffs(sigdists, tiles)
        cutoffs.update(tile_cutoffs)
        bases.update(tile_bases)

    write_outputs(cutoffs, bases)


# ===============
if is_snakemake_child:
    conf = params['conf']
    MIN_SPOTS_FOR_DIST = conf['polyA_seeder']['minimum_spots_for_dist_sampling']
    SPOTS_PREFER_NONLOCAL_SAMPLING = conf['polyA_seeder']['favored_spots_for_dist_sampling']
    KDE_BANDWIDTH = conf['polyA_seeder']['kde_bandwidth_for_dist']
    CUTOFF_RANGE_LOW = conf['polyA_seeder']['cutoff_score_search_low']
    CUTOFF_RANGE_HIGH = conf['polyA_seeder']['cutoff_score_search_high']
    TOTAL_CYCLES = params.total_cycles
    VERY_SMALL_PROBABILITY = 0.001
    NUM_THREADS = threads

    main()
