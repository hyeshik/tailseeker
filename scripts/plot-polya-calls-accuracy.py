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

from tailseeker.plotutils import colormap_lch
from tailseeker import stats
from functools import partial
import pandas as pd
import numpy as np
import os
from scipy.stats import mode

import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt

xtransform = lambda v: np.log2(v.clip(1)) if not isinstance(v, int) else np.log2(max(1, v))
xtransform_rev = lambda v: 2 ** v


def load_stats(options, controlsamples):
    res_stats = {}

    patbl = pd.read_csv(options.inputfile, index_col=0)

    dists = (patbl.cumsum(axis=0) / patbl.sum(axis=0))[[c for c, s in controlsamples]]
    dists.index = xtransform(np.array(dists.index))

    for sample, explength in controlsamples:
        polyadist = patbl[sample]
        polyadist_clipped = patbl[sample].iloc[1:].copy()
        polyadist_clipped.iloc[0] += polyadist.iloc[0]

        res_stats[sample] = {
            'median': stats.weighted_median(polyadist),
            'arimean': stats.weighted_mean(polyadist),
            'geomean': stats.weighted_geomean(polyadist_clipped),
            'mode': stats.weighted_mode(polyadist),
            'rmse': stats.weighted_rmse(polyadist, explength),
            'mae': stats.weighted_mae(polyadist, explength),
        }

    return dists, res_stats


def plot_dists(outpath, dists, controlsamples):
    from matplotlib import style

    style.use('ggplot')

    samplecolors = colormap_lch(len(controlsamples))
    samplessorted = sorted(controlsamples, key=lambda x: (x[1], x[0]))

    xticks = np.arange(0, dists.index[-1] + 0.1)
    xticks_disp = [format(v, 'g') for v in xtransform_rev(xticks)]

    fig, ax = plt.subplots(1, 1, figsize=(4, 3.2))
    ax.patch.set_facecolor('#f7f7f7')

    for (samplename, explength), color in zip(samplessorted, samplecolors):
        ax.plot(dists[samplename], color=color, label=samplename, zorder=3)
        ax.axvline(xtransform(explength), color=color, linewidth=2, alpha=0.3)

    leg = ax.legend(loc='center left', fontsize=10)
    plt.setp([leg.get_frame()], facecolor='white', edgecolor='#e8e8e8')

    ax.set_xlabel('Poly(A) length (nt)')
    ax.set_ylabel('Cumulative fraction')
    ax.set_title('Poly(A) length distribution', fontsize=12)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks_disp)

    plt.setp(ax.get_xgridlines() + ax.get_ygridlines(), color='#e0e0e0')

    plt.tight_layout()

    plt.savefig(outpath)

    plt.close(fig)


def write_descriptive_stats(outfile, stats, controlsamples):
    results = []

    for filepath, explength in controlsamples:
        samplename = os.path.basename(filepath).rsplit('.')[0]
        s = stats[filepath]

        results.append([
            samplename, int(explength), s['median'], s['mode'], s['geomean'],
            s['arimean'], s['rmse'], s['mae']])

    tbl = pd.DataFrame(results, columns=['name', 'design length', 'median',
                        'mode', 'geomean', 'arimean', 'rmse', 'mae'])

    for name, dtype in tbl.dtypes.items():
        if dtype == np.float64:
            tbl[name] = tbl[name].round(2)

    tbl.sort_values(by='design length').to_csv(outfile, index=False)

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description=
                            'Draw poly(A) length accuracy plots for poly(A) controls')
    parser.add_argument('--control', dest='control', metavar='LEN:NAME', type=str,
                        nargs='+')
    parser.add_argument('--input', dest='inputfile', metavar='FILE', type=str,
                        required=True)
    parser.add_argument('--output-plot', dest='output_plot', metavar='FILE', type=str,
                        default=None, help='Path to a PDF file')
    parser.add_argument('--output-stats', dest='output_stats', metavar='FILE', type=str,
                        default=None, help='Path to a CSV file')

    options = parser.parse_args()

    controlsamples = [
        (name, int(length))
        for length, name in map(partial(str.split, sep=':'), options.control)
    ]

    return options, controlsamples


if __name__ == '__main__':
    import pickle

    options, controlsamples = parse_arguments()

    dists, stats = load_stats(options, controlsamples)

    if options.output_plot is not None:
        plot_dists(options.output_plot, dists, controlsamples)

    if options.output_stats is not None:
        write_descriptive_stats(options.output_stats, stats, controlsamples)

