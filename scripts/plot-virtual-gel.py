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

from tailseeker.stats import gaussian_kde
from itertools import chain
import pandas as pd
import numpy as np

import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import AutoMinorLocator


# TODO: generate list of marker positions in adaptive way to the range.
MARKER_POSITIONS = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 235]

LANE_WIDTH = 7
LANE_SPACE = 1
LANE_BODY = LANE_WIDTH - LANE_SPACE
LANE_CONTROLGUIDE = 0.15

IMAGE_HEIGHT = 200



def construct_virtual_gel_image(tagdensities, controlsamples, expsamples, height=IMAGE_HEIGHT):
    tagdensities = tagdensities[controlsamples + expsamples].copy()
    #tagdensities = tagdensities[expsamples + controlsamples].copy()

    is_control = tagdensities.columns.map(lambda x: x in controlsamples)

    if is_control.sum() > 0:
        # Rescale the control samples to the maximum in each lane.
        control_signals = tagdensities.loc[:, is_control]
        tagdensities.loc[:, is_control] = control_signals / control_signals.max(axis=0)

    # Rescale the control samples to the total poly(A)+ tag count in each lane.
    experimental_signals = tagdensities.loc[:, ~is_control]
    tagdensities.loc[:, ~is_control] = experimental_signals / experimental_signals.max().max()

    # Pad the gel image with signals and spacings between lanes.
    virtual_width = tagdensities.shape[1] * LANE_WIDTH + LANE_SPACE
    imgdata = np.zeros((height, virtual_width), np.float64)

    pos = 0
    for i, (vpos, td) in enumerate(tagdensities.items()):
        if i == 0:
            pos += LANE_SPACE # vertical spacing on the left

        for j in range(LANE_BODY):
            imgdata[:, pos] = td
            pos += 1

        pos += LANE_SPACE

    return imgdata


def size2phys_fun(min_polya, max_polya):
    migrfun = lambda x: x ** 0.5 - 1
    phys_low = migrfun(min_polya)
    phys_width = migrfun(max_polya) - phys_low

    def size2phys(x):
        return (migrfun(x) - phys_low) / phys_width

    return size2phys


def plot(options, controlsamples, expsamples):
    # Load tag counts by poly(A) lengths
    tagcounts = pd.read_csv(options.tagcounts, index_col=0)
    tagcounts_pa = tagcounts.iloc[options.min_polya_len:]
    if options.normalize_by_total_tags:
        tagcounts_pa = tagcounts_pa.divide(tagcounts.sum(axis=0))
    else:
        tagcounts_pa = tagcounts_pa.divide(tagcounts_pa.sum(axis=0))

    if options.merge_controls and controlsamples:
        ctlsum = tagcounts_pa[controlsamples].sum(axis=1)
        newbc = tagcounts_pa[expsamples].copy()
        newbc['ctl'] = ctlsum
        tagcounts_pa = newbc
        controlsamples = ['ctl']

    samples = controlsamples + expsamples

    # Calculate kernel density estimates
    max_polya_len = int(np.ceil(tagcounts_pa.index.max() / 5) * 5)
    size2phys = size2phys_fun(options.min_polya_len, max_polya_len)
    vpositions = np.linspace(size2phys(options.min_polya_len),
                             size2phys(max_polya_len), IMAGE_HEIGHT)

    kdes = pd.DataFrame({
        name: gaussian_kde(
                list(size2phys(np.array(cnts.index))),
                weights=list(cnts),
                bw_method=options.kde_bandwidth)(vpositions)
        for name, cnts in tagcounts_pa.items()}, index=vpositions)

    # Construct the image and plot it.
    imgdata = construct_virtual_gel_image(kdes, controlsamples, expsamples)

    figwidth = options.width if options.width is not None else (
                                       1.4 + 0.33 * len(samples))
    fig, ax = plt.subplots(1, 1, figsize=(figwidth, options.height))

    ax.pcolor(imgdata, cmap='Greys', rasterized=True, vmin=0, vmax=1)

    sample_label_tickpos = np.arange(
        LANE_SPACE + (LANE_WIDTH - LANE_SPACE) / 2,
        LANE_SPACE + LANE_WIDTH * len(samples), LANE_WIDTH)

    ax.set_xticks(sample_label_tickpos)
    ax.set_xticklabels(samples, rotation=60, ha='left', va='bottom')

    ax.xaxis.tick_top()

    ax.set_ylabel('Poly(A) tail length (nt)')
    ax.set_yticks([(size2phys(m) - size2phys(options.min_polya_len)) /
                   (size2phys(max_polya_len) - size2phys(options.min_polya_len))
                    * IMAGE_HEIGHT for m in MARKER_POSITIONS])
    ax.set_yticklabels(MARKER_POSITIONS)

    for spinepos in 'top bottom left right'.split():
        ax.spines[spinepos].set_visible(False)

    plt.setp(ax.get_xticklines(), visible=False)
    ax.tick_params(labelright=True)
    ax.tick_params(axis='y', direction='out')

    ax.set_xlim(0, imgdata.shape[1])
    ax.set_ylim(size2phys(options.min_polya_len) * IMAGE_HEIGHT, imgdata.shape[0])

    plt.tight_layout()

    plt.savefig(options.output_plot, dpi=200)


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description=
                            'Draw a virtual gel image for poly(A) length distributiions')
    parser.add_argument('--tagcounts', dest='tagcounts', metavar='FILE', type=str,
                        required=True,
                        help='A CSV file containing tag counts for poly(A) lengths')
    parser.add_argument('--controls', dest='controls', metavar='NAME,NAME..', type=str,
                        help='Names of poly(A) length control samples to show '
                             '(comma-separated)')
    parser.add_argument('--samples', dest='samples', metavar='NAME,NAME..', type=str,
                        help='Names of experimental samples to show (comma-separated)')
    parser.add_argument('--width', dest='width', metavar='INCHES', type=float,
                        help='Width of the full image (in inches)')
    parser.add_argument('--height', dest='height', metavar='INCHES', type=float,
                        default=6., help='Height of the full image (in inches)')
    parser.add_argument('--output-plot', dest='output_plot', metavar='FILE', type=str,
                        required=True, help='Output path of the resulting gel image')
    parser.add_argument('--by-total-tag-counts', dest='normalize_by_total_tags',
                        action='store_true', default=False,
                        help='Normalize tag counts by total (default: by poly(A)+ tags)')
    parser.add_argument('--minimum-polya-length', dest='min_polya_len', type=int,
                        default=5, help='Poly(A) length at the bottom of the gel')
    parser.add_argument('--kde-bandwidth', dest='kde_bandwidth',
                        metavar='VALUE', type=float, default=0.07,
                        help='Bandwidth for the kernel density estimation')
    parser.add_argument('--merge-controls', dest='merge_controls',
                        action='store_true', default=False,
                        help='Merge control lanes into one to form a ladder')

    options = parser.parse_args()

    controlsamples = options.controls.split(',') if options.controls else []
    expsamples = options.samples.split(',') if options.samples else []

    return options, controlsamples, expsamples


if __name__ == '__main__':
    import pickle

    options, controlsamples, expsamples = parse_arguments()

    plot(options, controlsamples, expsamples)

