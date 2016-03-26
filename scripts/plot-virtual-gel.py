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

from itertools import chain
import pandas as pd
import numpy as np

import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import AutoMinorLocator


POLYA_BIN_BREAKS = [
    1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 28, 31, 34, 37, 40,
    43, 46, 49, 52, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99,
    104, 109, 114, 119, 124, 129, 134, 139, 144, 149, 154, 160, 166,
    172, 178, 184, 190, 196, 202, 208, 214, 220, 226, 235,
]

MARKER_POSITIONS = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 235]
size2phys = lambda x: (x ** 0.5 - 1 ** 0.5) / (235 ** 0.5 - 1 ** 0.5)
MINIMUM_POLYA_LENGTH = 3

IMAGE_HEIGHT = 1000
SMOOTHEDGE = [.25, .5, .85]
LANE_WIDTH = 31
LANE_SMOOTHEDGE = len(SMOOTHEDGE)
LANE_SPACE = 2
LANE_BODY = LANE_WIDTH - LANE_SMOOTHEDGE * 2 - LANE_SPACE
LANE_CONTROLGUIDE = 0.15


def construct_virtual_gel_image(tagcounts, controlsamples, expsamples, height=IMAGE_HEIGHT):
    tagcounts = tagcounts[controlsamples + expsamples]

    # Adjust intensities to make tag populations proportional to area*intensity.
    adjusted_values = []
    covering_ranges = []

    for idx, row in tagcounts.iterrows():
        low = int(size2phys(POLYA_BIN_BREAKS[idx]) * height)
        high = int(size2phys(POLYA_BIN_BREAKS[idx+1]) * height)
        density = row / (high - low)

        covering_ranges.append((low, high))
        adjusted_values.append(density)

    adjusted_values = np.array(adjusted_values)

    is_control = tagcounts.columns.map(lambda x: x in controlsamples)

    # Rescale the control samples to the maximum in each lane.
    control_signals = adjusted_values[:, is_control]
    adjusted_values[:, is_control] = control_signals / control_signals.max(axis=0)

    # Rescale the control samples to the total poly(A)+ tag count in each lane.
    experimental_signals = adjusted_values[:, ~is_control]
    adjusted_values[:, ~is_control] = experimental_signals / experimental_signals.max().max()

    # Pad the gel image with signals and spacings between lanes.
    virtual_width = tagcounts.shape[1] * LANE_WIDTH + LANE_SPACE
    imgdata = np.zeros((height, virtual_width), np.float64)

    signalmask = np.array(SMOOTHEDGE + [1] * LANE_BODY + SMOOTHEDGE[::-1] +
                          [0] * LANE_SPACE)
    signalboost = [0] * LANE_WIDTH
    signalboost[LANE_WIDTH // 2] = LANE_CONTROLGUIDE

    for idx, signals, (low, high) in zip(tagcounts.index, adjusted_values, covering_ranges):
        densityrow = list(chain(*([[0.] * LANE_SPACE] +
                                  [list(np.max([v * signalmask, signalboost], axis=0))
                                   if is_ctl else list(v * signalmask)
                                   for is_ctl, v in zip(is_control, signals)])))
        imgdata[low:high] = densityrow

    return imgdata


def plot(options, controlsamples, expsamples):
    tagcounts = pd.read_csv(options.tagcounts, index_col=0)
    tagcounts_pa = tagcounts.iloc[MINIMUM_POLYA_LENGTH:]
    if options.normalize_by_total_tags:
        tagcounts_pa = tagcounts_pa.divide(tagcounts.sum(axis=0))
    else:
        tagcounts_pa = tagcounts_pa.divide(tagcounts_pa.sum(axis=0))

    # Aggregate the tag counts in the fake gel bins.
    binnedcounts = tagcounts_pa.copy()
    binnedcounts['bin'] = pd.cut(tagcounts_pa.index, POLYA_BIN_BREAKS, right=False,
                                 include_lowest=True, labels=False)
    binnedcounts = binnedcounts.groupby('bin').agg(np.sum)
    binnedcounts.index = list(map(int, binnedcounts.index))

    samples = controlsamples + expsamples


    # Construct the image and plot it.
    imgdata = construct_virtual_gel_image(binnedcounts, controlsamples, expsamples)

    figwidth = options.width if options.width is not None else (
                                       1.4 + 0.33 * len(samples))

    fig = plt.figure(figsize=(figwidth, options.height))
    ax = plt.gca()

    ax.pcolor(imgdata, cmap='Greys', rasterized=True, vmin=0, vmax=1)

    sample_label_tickpos = np.arange(
        LANE_SPACE + (LANE_WIDTH - LANE_SPACE) / 2,
        LANE_SPACE + LANE_WIDTH * len(samples), LANE_WIDTH)

    ax.set_xticks(sample_label_tickpos)
    ax.set_xticklabels(samples, rotation=60, ha='left', va='bottom')

    ax.xaxis.tick_top()

    ax.set_ylabel('Poly(A) tail length (nt)')
    ax.set_yticks([int(size2phys(m) * IMAGE_HEIGHT) for m in MARKER_POSITIONS])
    ax.set_yticklabels(MARKER_POSITIONS)

    for spinepos in 'top bottom left right'.split():
        ax.spines[spinepos].set_visible(False)

    plt.setp(ax.get_xticklines(), visible=False)
    ax.tick_params(labelright=True)
    ax.tick_params(axis='y', direction='out')

    ax.set_xlim(0, imgdata.shape[1])
    ax.set_ylim(size2phys(MINIMUM_POLYA_LENGTH) * IMAGE_HEIGHT, imgdata.shape[0])

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

    options = parser.parse_args()

    controlsamples = options.controls.split(',') if options.controls else []
    expsamples = options.samples.split(',') if options.samples else []

    return options, controlsamples, expsamples


if __name__ == '__main__':
    import pickle

    options, controlsamples, expsamples = parse_arguments()

    plot(options, controlsamples, expsamples)

