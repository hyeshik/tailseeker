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

from tailor.parsers import parse_polya_calls
from tailor.plotutils import prepare_cumulative, colormap_lch
from functools import partial
import pandas as pd
import numpy as np
import os
from scipy.stats import mode

import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import AutoMinorLocator

xtransform = lambda v: np.log2(v.clip(1)) if not isinstance(v, int) else np.log2(max(1, v))
xtransform_rev = lambda v: 2 ** v

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


def geomean(lengths):
    return np.exp(np.log(lengths.clip(1)).mean())

def rmsd(lengths, explength):
    return ((lengths - explength) ** 2).mean() ** 0.5

def load_stats(options, controlsamples):
    stats = {}
    dists = {}

    for filepath, explength in controlsamples:
        patbl = parse_polya_calls.as_table(filepath)

        finaldist_x, finaldist_y = prepare_cumulative(xtransform(patbl['polya_len']))
        seqdist_x, seqdist_y = prepare_cumulative(xtransform(patbl['seqbased_len']))
        hmmdist_x, hmmdist_y = prepare_cumulative(xtransform(patbl['hmmbased_len']))
        naivedist_x, naivedist_y = prepare_cumulative(xtransform(patbl['naive_len']))

        dists[filepath] = {
            'final': (finaldist_x, finaldist_y),
            'seqbased': (seqdist_x, seqdist_y),
            'hmmbased': (hmmdist_x, hmmdist_y),
            'naive': (naivedist_x, naivedist_y),
        }

        stats[filepath] = {
            'median_final': np.median(patbl['polya_len']),
            'median_seqbased': np.median(patbl['seqbased_len']),
            'median_hmmbased': np.median(patbl['hmmbased_len']),
            'median_naive': np.median(patbl['naive_len']),

            'arimean_final': np.mean(patbl['polya_len']),
            'arimean_seqbased': np.mean(patbl['seqbased_len']),
            'arimean_hmmbased': np.mean(patbl['hmmbased_len']),
            'arimean_naive': np.mean(patbl['naive_len']),

            'geomean_final': geomean(patbl['polya_len']),
            'geomean_seqbased': geomean(patbl['seqbased_len']),
            'geomean_hmmbased': geomean(patbl['hmmbased_len']),
            'geomean_naive': geomean(patbl['naive_len']),

            'mode_final': mode(patbl['polya_len']).mode[0],
            'mode_seqbased': mode(patbl['seqbased_len']).mode[0],
            'mode_hmmbased': mode(patbl['hmmbased_len']).mode[0],
            'mode_naive': mode(patbl['naive_len']).mode[0],

            'rmsd_final': rmsd(patbl['polya_len'], explength),
            'rmsd_seqbased': rmsd(patbl['seqbased_len'], explength),
            'rmsd_hmmbased': rmsd(patbl['hmmbased_len'], explength),
            'rmsd_naive': rmsd(patbl['naive_len'], explength),
        }

    return dists, stats


def plot_dists(outdir, dists, controlsamples):
    from matplotlib import style

    style.use('ggplot')

    samplecolors = colormap_lch(len(controlsamples))
    samplessorted = sorted(controlsamples, key=lambda x: (x[1], x[0]))

    xtickmax = np.ceil(max(x[-1] for d1 in dists.values() for x, y in d1.values()))
    xticks = np.arange(0, xtickmax + 0.1)
    xticks_disp = [format(v, 'g') for v in xtransform_rev(xticks)]

    for plottype in next(iter(dists.values())).keys():
        fig = plt.figure(figsize=(4, 3.2))

        ax = fig.add_subplot(1, 1, 1)

        for (filepath, explength), color in zip(samplessorted, samplecolors):
            samplename = os.path.basename(filepath).rsplit('.')[0]
            distx, disty = dists[filepath][plottype]
            ax.plot(distx, disty, color=color, label=samplename, zorder=3)
            ax.axvline(xtransform(explength), color=color, linewidth=2, alpha=0.3)

        ax.legend(loc='center left', fontsize=10)
        ax.set_xlabel('Poly(A) length (nt)')
        ax.set_ylabel('Cumulative fraction')
        ax.set_title('Poly(A) control distribution: {}'.format(plottype), fontsize=12)

        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_disp)

        plt.tight_layout()

        output_path = os.path.join(outdir, 'control-length-accuracy-{}.pdf'.format(plottype))
        plt.savefig(output_path)

        plt.close(fig)


def write_descriptive_stats(outfile, stats, controlsamples):
    results = []

    for filepath, explength in controlsamples:
        samplename = os.path.basename(filepath).rsplit('.')[0]
        s = stats[filepath]

        results.append([
            samplename, explength,
            s['median_final'], s['mode_final'], s['geomean_final'],
            s['arimean_final'], s['rmsd_final'],
            s['median_seqbased'], s['mode_seqbased'],
            s['geomean_seqbased'], s['arimean_seqbased'], s['rmsd_seqbased'],
            s['median_hmmbased'], s['mode_hmmbased'],
            s['geomean_hmmbased'], s['arimean_hmmbased'], s['rmsd_hmmbased'],
            s['median_naive'], s['mode_naive'], s['geomean_naive'], s['arimean_naive'],
            s['rmsd_naive'],
        ])

    tbl = pd.DataFrame(results, columns=['name', 'design length', 'median (final)',
                        'mode (final)', 'geomean (final)', 'arimean (final)',
                        'rmsd (final)', 'median (seq)', 'mode (seq)', 'geomean (seq)',
                        'arimean (seq)', 'rmsd (seq)', 'median (HMM)', 'mode (HMM)',
                        'geomean (HMM)', 'arimean (HMM)', 'rmsd (HMM)', 'median (naive)',
                        'mode (naive)', 'geomean (naive)', 'arimean (naive)',
                        'rmsd (naive)'])

    for name, dtype in tbl.dtypes.items():
        if dtype == np.float64:
            tbl[name] = tbl[name].round(2)

    tbl.to_csv(outfile, index=False)

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description=
                            'Draw poly(A) length accuracy plots for poly(A) controls')
    parser.add_argument('--control', dest='control', metavar='LEN:FILE', type=str,
                        nargs='+')
    parser.add_argument('--output-plots', dest='output_plots', metavar='DIR', type=str,
                        default=None, help='Directory to write plot files')
    parser.add_argument('--output-stats', dest='output_stats', metavar='FILE', type=str,
                        default=None, help='Path to a CSV file')

    options = parser.parse_args()

    controlsamples = [
        (filepath, int(length))
        for length, filepath
        in map(partial(str.split, sep=':'), options.control)
    ]

    return options, controlsamples


if __name__ == '__main__':
    import pickle

    options, controlsamples = parse_arguments()

    dists, stats = load_stats(options, controlsamples)

    if options.output_plots is not None:
        plot_dists(options.output_plots, dists, controlsamples)

    if options.output_stats is not None:
        write_descriptive_stats(options.output_stats, stats, controlsamples)

