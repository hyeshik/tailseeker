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

from tailseeker.plotutils import apply_dropped_spine, prepare_cumulative
from tailseeker import stats
from struct import unpack
import pandas as pd
import numpy as np
import random
from scipy import interpolate
import gzip
from glob import glob
import sys
import re

import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import cm


def load_sigdump(filename):
    inpfile = gzip.open(filename, 'rb')
    ncycles, datasize = unpack('<II', inpfile.read(8))
    sigs = np.fromstring(inpfile.read(), np.float32)
    return sigs.reshape((sigs.shape[0] // ncycles, ncycles))

def load_all_sigdumps(sigdumpfiles):
    files = glob(sigdumpfiles.replace('%%TILE%%', '*'))
    pathpattern = re.compile(sigdumpfiles.replace('%%TILE%%', '([^.-]+)'))
    loaded_dumps = {}

    print('Loading signal dumps... ', end='', file=sys.stderr)
    sys.stderr.flush()

    for filepath in sorted(files):
        tilenumber = pathpattern.findall(filepath)
        if not tilenumber:
            continue

        loaded_dumps[tilenumber[0]] = load_sigdump(filepath)

    print('{num_spots} spots from {num_tiles} tiles '
          '({min_spots_per_tile}-{max_spots_per_tile} spots per tile)'.format(
            num_spots=sum(map(len, loaded_dumps.values())),
            num_tiles=len(loaded_dumps),
            min_spots_per_tile=min(map(len, loaded_dumps.values())),
            max_spots_per_tile=max(map(len, loaded_dumps.values()))), file=sys.stderr)

    return loaded_dumps

def plot_signal_hist(sigdumps, options):
    print('Drawing a histogram for signal scores...', file=sys.stderr)

    binfactor = 100
    saturation_percentile = 98
    bins = binfactor * options.signal_max

    r = []
    for i in range(sigdumps.shape[1]):
        r.append(np.histogram(sigdumps[:, i], range=(0, options.signal_max),
                              bins=bins)[0])
    histp = np.array(r)
    vmax = np.percentile(histp.ravel(), saturation_percentile)

    fig, ax = plt.subplots(1, 1, figsize=(7, 3.5))

    ax.pcolor(histp.T, vmax=vmax, cmap='Blues')
    ax.axvline(options.design_length, c='red')

    ax.set_xlabel('Poly(A) cycle')
    ax.set_ylabel('Score')
    ax.set_xlim(0, sigdumps.shape[1])
    ax.set_ylim(0, 145)
    ticks = np.arange(0, options.signal_max+0.001, 0.2)
    ax.set_yticks(ticks * binfactor)
    ax.set_yticklabels(ticks)

    apply_dropped_spine(ax, xgrid=True)

    plt.tight_layout()
    plt.savefig(options.signal_hist_output)
    plt.close(fig)

def determine_threshold(sigdumps, options):
    print('Finding the optimal score thresholds for poly(A)... ', end='', file=sys.stderr)
    sys.stderr.flush()

    polya_region = slice(0, options.design_length - options.polya_spacing + 1)
    genebody_region = slice(options.design_length + options.genebody_spacing,
                            sigdumps.shape[1])
    scoring_bins = 2500
    pseudocount = 1

    cnt_pa, edges_pa = np.histogram(sigdumps[:, polya_region].ravel(),
                                    range=(0, options.signal_max), bins=scoring_bins)
    cnt_nonpa, edges_nonpa = np.histogram(sigdumps[:, genebody_region].ravel(),
                                    range=(0, options.signal_max), bins=scoring_bins)

    pa_pdf = (cnt_pa + pseudocount) / cnt_pa.sum()
    nonpa_pdf = (cnt_nonpa + pseudocount) / cnt_nonpa.sum()
    logLR = stats.smooth(np.log(pa_pdf / nonpa_pdf), options.pdf_smoothing)

    bincenters = (edges_pa[:-1] + edges_pa[1:]) / 2

    if options.preset_threshold is None:
        expected_range = bincenters < options.expected_maximum

        inwin_x = bincenters[expected_range]
        inwin_y = logLR[expected_range]
        inwin_f = interpolate.UnivariateSpline(inwin_x, inwin_y, s=0)

        cutoff_val = inwin_f.roots()[0]
    else:
        cutoff_val = options.preset_threshold

    if options.param_det_output is not None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 3.5))

        ax.plot(bincenters, logLR, c='#1040a0')
        ax.axvline(cutoff_val, c='red')
        ax.axhline(0, c='gray')
        ax.annotate('{:.4f}'.format(cutoff_val), (cutoff_val, 0),
                    xytext=(5, -5), textcoords='offset points', va='top')

        ax.set_xlim(0, options.expected_maximum)
        ax.set_xlabel('Score')
        ax.set_ylabel('Log likelihood ratio')

        apply_dropped_spine(ax, xgrid=True)

        plt.tight_layout()
        plt.savefig(options.param_det_output)

        plt.close(fig)

    print('->', format(cutoff_val, '.4f'), file=sys.stderr)

    return cutoff_val

def try_measurements(sigdumps, cutoff_val, options):
    print('Measuring poly(A) tail lengths in draft quality...', file=sys.stderr)
    draft_palens = ((sigdumps > cutoff_val) * 2 - 1).cumsum(axis=1).argmax(axis=1)

    # Apply downhill-extension
    choose = lambda a, c: np.array([c[a[I]][I] for I in np.ndindex(a.shape)])

    maxindex = sigdumps.shape[1] - 1
    downhillarray = np.array([choose((draft_palens + i).clip(0, maxindex), sigdumps.T)
                              for i in range(-1, 10)]).T
    zeropad = np.array([0] * len(draft_palens))[np.newaxis].T
    slopedown = np.append((np.diff(downhillarray) < 0).astype('i'), zeropad, axis=1)
    downhill_ext = np.ceil(slopedown.argmin(axis=1) / 2).astype('i')

    palens = draft_palens + downhill_ext

    if options.signal_samples_output is not None:
        print("Writing a plot for signal samples with the optimized threshold applied...",
              file=sys.stderr)
        random_samples = 3000
        sampled_indices = list(range(len(sigdumps)))
        random.shuffle(sampled_indices)
        del sampled_indices[random_samples:]

        samples = sigdumps[sampled_indices]
        samples_order = np.argsort(-palens[sampled_indices])

        fig, ax = plt.subplots(1, 1, figsize=(7, 5.0))

        pseudoscore = 0.00001
        samples_transformed = np.log2(
            (samples[samples_order] + pseudoscore) / (cutoff_val + pseudoscore))
        ax.pcolor(samples_transformed, cmap='PiYG_r', vmax=3, vmin=-3, rasterized=True)
        ax.axvline(options.design_length, c='#a03030', alpha=.6)
        ax.set_xlim(0, sigdumps.shape[1])
        ax.set_xlabel('Poly(A) cycle')
        ax.set_ylabel('Spot samples')
        ax.set_yticks(np.linspace(0, random_samples, 11))
        ax.set_yticklabels(np.arange(0, 101, 10))
        ax.invert_yaxis()

        plt.setp(ax.get_yticklines(), visible=False)

        plt.tight_layout()
        plt.savefig(options.signal_samples_output)

        plt.close(fig)

    if options.cdf_plot_output is not None:
        print("Writing a CDF plot for poly(A) measurements...", file=sys.stderr)

        cdf_x, cdf_y = prepare_cumulative(palens)

        fig, ax = plt.subplots(1, 1, figsize=(4.5, 3.5))

        ax.plot(cdf_x, cdf_y, c='#3030a0')
        ax.axvline(options.design_length, c='red')
        ax.set_xlabel('Poly(A) cycle')
        ax.set_ylabel('Cumulative fraction')

        apply_dropped_spine(ax, xgrid=True)

        plt.tight_layout()
        plt.savefig(options.cdf_plot_output)

        plt.close(fig)

    if options.measurement_csv_output is not None:
        print("Writing a count table for draft measurements...", file=sys.stderr)

        vcounts = pd.Series(palens).value_counts()
        pd.DataFrame({'spots': vcounts}).sort_index().to_csv(options.measurement_csv_output)


def main(options):
    sigdumps_bytile = load_all_sigdumps(options.sigdumpfiles)
    sigdumps = np.vstack(sigdumps_bytile.values())

    if options.signal_hist_output is not None:
        plot_signal_hist(sigdumps, options)

    cutoff_val = determine_threshold(sigdumps, options)

    try_measurements(sigdumps, cutoff_val, options)

    print('{:.6f}'.format(cutoff_val))


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description=
                            'Draw a virtual gel image for poly(A) length distributiions')
    parser.add_argument('--sigdump-files', dest='sigdumpfiles', type=str,
                        required=True, help='Path to sigdump files with a %%TILE%% placeholder')
    parser.add_argument('--output-signal-histogram', dest='signal_hist_output',
                        type=str, default=None,
                        help='Path to write a plot for signal histogram')
    parser.add_argument('--output-parameter-determination', dest='param_det_output',
                        type=str, default=None,
                        help='Path to write a plot for the parameter determination')
    parser.add_argument('--output-signal-samples', dest='signal_samples_output',
                        type=str, default=None,
                        help='Path to write a plot for the binarized signal samples')
    parser.add_argument('--output-measurements-cdf-plot', dest='cdf_plot_output',
                        type=str, default=None,
                        help='Path to write a plot for a poly(A) measurement CDF.')
    parser.add_argument('--output-measurements-csv', dest='measurement_csv_output',
                        type=str, default=None,
                        help='Path to write a CSV table for a poly(A) measurement CDF.')
    parser.add_argument('--signal-max', dest='signal_max', type=float,
                        default=1.45, help='Maximum signal value')
    parser.add_argument('--design-length', dest='design_length', type=float,
                        default=118, help='Length of poly(A) tail in design')
    parser.add_argument('--polya-region-stepback', dest='polya_spacing', type=int,
                        default=5, help='Distance from the design length to the end of '
                                        'poly(A) region to be sampled')
    parser.add_argument('--genebody-region-stepforward', dest='genebody_spacing', type=int,
                        default=5, help='Distance from the design length to the begin of '
                                        'gene body region to be sampled')
    parser.add_argument('--pdf-smoothing-window', dest='pdf_smoothing', type=int,
                        default=5, help='Window size for smoothing PDFs')
    parser.add_argument('--expected-maximum', dest='expected_maximum', type=float,
                        default=0.2, help='Empirically expected maximum value for a threshold')
    parser.add_argument('--preset-threshold', dest='preset_threshold', type=float,
                        default=None, help='Evaluate a given parameter rather than '
                                           'producing an optimal value')

    return parser.parse_args()


if __name__ == '__main__':
    options = parse_arguments()
    main(options)

