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

from tailor.parsers import parse_pascore, parse_sqi
from tailor.stats import sample_iterable, similarity_sort
from tailor.fileutils import ParallelMatchingReader, open_gzip_buffered
from rpy2.robjects import numpy2ri; numpy2ri.activate()
from rpy2 import robjects
import numpy as np
from Bio.Cluster import cluster
import random
import sys
import subprocess as sp
from operator import itemgetter

import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import AutoMinorLocator


def load_pascore_trimmed(pascorefile, trim_length):
    for spot in parse_pascore(pascorefile):
        pascore = spot.pascore
        if len(pascore) >= trim_length:
            yield '{}:{:08d}'.format(spot.tile, spot.cluster), pascore[:trim_length]


def run_outlier_removal(options):
    r = robjects.r
    r.library('mvoutlier', quietly=True)
    r.library('robustbase', quietly=True)
    r("""\
    pickoutliers <- function (x, quan={options.quan}, alpha={options.alpha}) {{
        rob <- covMcd(x, alpha=quan)
        xarw <- arw(x, rob$center, rob$cov, alpha=alpha)
        dist <- mahalanobis(x, center=rob$center, cov=rob$cov)
        o <- (sqrt(dist) > min(sqrt(xarw$cn), sqrt(qchisq(1-alpha, dim(x)[2]))))
        l <- list(outliers=o, md=sqrt(dist))
        l
    }}""".format(options=options))


    # random sample for pass 1
    parserfeed = load_pascore_trimmed(options.input_pascore, options.trim)

    sampled_lines = sample_iterable(parserfeed, options.pass1)
    trimlength = min([len(sl) for name, sl in sampled_lines])
    if options.trim is not None and options.trim < trimlength:
        trimlength = options.trim
    pass1names = [name for name, row in sampled_lines]
    pass1spots = np.array([row[:trimlength] for name, row in sampled_lines])

    # trim out preamble and uneven right end for granule size
    if pass1spots.shape[1] % options.granule > 0:
        p1trimmed = pass1spots[:, :pass1spots.shape[1] // options.granule * options.granule]
    else:
        p1trimmed = pass1spots

    # merge adjacent cells to stabilize the signal for outlier detection
    p1merged = p1trimmed.reshape([p1trimmed.shape[0], p1trimmed.shape[1] // options.granule,
                                  options.granule]).sum(axis=2)

    # check outliers and remove outliers
    outliers = np.array(list(r.pickoutliers(p1merged).rx('outliers')[0]))
    nonoutliers = np.array(1 - outliers, 'bool')

    survivors = [(name, spot)
                 for name, spot, insider in zip(pass1names, pass1spots, nonoutliers)
                 if insider]

    # pass 2 - get random samples from the survivors from pass 1
    if len(survivors) < options.pass2:
        print("Not enough samples for pass 3 ({} < {})".format(
                                len(survivors), options.pass2), file=sys.stderr)
        sys.exit(1)

        assert len(survivors) >= options.pass2
    random.shuffle(survivors)
    del survivors[options.pass2:]

    survivors.sort(key=itemgetter(0))

    return np.array(pass1spots), np.array([s[1] for s in survivors]), \
           [s[0] for s in survivors]


def run_specific_picking(options):
    readpool = open(options.readids).read().split()
    random.shuffle(readpool)
    del readpool[options.pass1:]

    scores = []
    for readid, pascore in load_pascore_trimmed(options.input_sqi, options.input_pascore, withid=True):
        if readid in readpool:
            scores.append(pascore)
            if len(scores) >= options.pass2:
                break

    trimlength = min(map(len, scores))
    if options.trim is not None and options.trim < trimlength:
        trimlength = options.trim

    finalspots = np.array([row[:trimlength] for row in scores])

    return finalspots, list(finalspots)


def plot_qc(p1result, p2result, output, clip_range):
    fig = plt.figure(figsize=(8, 4.5))
    clip_range = list(map(float, clip_range.split(':')))

    def plot_an_array(arr):
        plt.pcolor(np.array(arr).clip(*clip_range), cmap=cm.OrRd, rasterized=True,
                    vmin=clip_range[0], vmax=clip_range[1])
        #plt.axvline(25, c='black')
        #plt.axvline(128, c='black')
        plt.grid(True)
        plt.grid(True, which='minor')
        plt.gca().yaxis.grid(False)
        plt.gca().yaxis.grid(False, which='minor')
        #plt.colorbar()
        plt.ylim(ymax=len(arr))
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
        plt.xlabel('cycle (read 2)')
        plt.ylabel('reads')

    fig.add_subplot(1, 2, 1)
    random.shuffle(p1result)
    neworder = similarity_sort(p1result[:len(p2result)])
    plot_an_array(np.array(p1result)[neworder])

    fig.add_subplot(1, 2, 2)
    neworder = similarity_sort(p2result)
    plot_an_array(np.array(p2result)[neworder])

    plt.tight_layout()

    plt.savefig(output)

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Pick random samples to train models '
                                                 'from spike-in clusters without outliers')
    parser.add_argument('--input-pascore', dest='input_pascore', metavar='FILE', type=str,
                        required=True, help='Path to a file containing PA scores')
    parser.add_argument('--input-sqi', dest='input_sqi', metavar='FILE', type=str,
                        default=None, help='Path to a SQI file')
    parser.add_argument('--pass1', dest='pass1', metavar='N', type=int,
                        help='Number of samples in the first pass (before outlier removal)',
                        default=2000)
    parser.add_argument('--pass2', dest='pass2', metavar='N', type=int,
                        help='Number of samples in the second pass (after outlier removal)',
                        default=1000)
    parser.add_argument('--quan', dest='quan', metavar='alpha', type=float,
                        help='Alpha value for covariance MCD', default=0.5)
    parser.add_argument('--alpha', dest='alpha', metavar='alpha', type=float,
                        help='Alpha value for arw', default=0.025)
    parser.add_argument('--granule-size', dest='granule', metavar='size', type=int,
                        help='Size of bins for signal averaging.', default=15)
    parser.add_argument('--read-ids', dest='readids', metavar='FILE', type=str,
                        default=None, help='Specify read IDs instead of outlier removal')
    parser.add_argument('--trim', dest='trim', metavar='N', type=int,
                        required=True, help='Number of cycles to keep from the first')
    parser.add_argument('--output', dest='output', metavar='FILE', type=str,
                        default=None, help='Path to a file containing PA scores')
    parser.add_argument('--output-qc-plot', dest='output_qc_plot', metavar='FILE', type=str,
                        default=None, help='Path to put a plot for QC')
    parser.add_argument('--qc-plot-range', dest='qc_plot_range', metavar='MIN:MAX', type=str,
                        default='-2:2', help='Value range for a QC plot')
    parser.add_argument('--output-training-set-list', dest='output_training_set_list',
                        metavar='FILE', type=str, default=None,
                        help='Path to write the list of read ids included in training set')

    options = parser.parse_args()

    return options


if __name__ == '__main__':
    import pickle

    options = parse_arguments()
    if options.readids is None:
        p1result, p2result, p2names = run_outlier_removal(options)
    else:
        p1result, p2result = run_specific_picking(options)

    if options.output is not None:
        np.save(options.output, p2result)

    if options.output_qc_plot is not None:
        plot_qc(p1result, p2result, options.output_qc_plot, options.qc_plot_range)

    if options.output_training_set_list and options.readids is None:
        output = open(options.output_training_set_list, 'w')
        print('\n'.join(p2names), file=output)

