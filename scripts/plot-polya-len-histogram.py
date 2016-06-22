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
from tailseeker.plotutils import apply_dropped_spine, colormap_lch
import pandas as pd
import numpy as np

import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt

MINIMUM_POLYA_LEN = snakemake.params.minimum_polya_length
KDE_BANDWIDTH = snakemake.params.kde_bandwidth
XSCALE_TRANSFORM_FACTOR = snakemake.params.x_transform_factor
EXCLUDE_SAMPLES = snakemake.params.exclude

xticks = pd.Series([5] + list(np.arange(10, 100, 10)) + list(np.arange(100, 500, 25)))

polyadist = pd.read_csv(snakemake.input[0], index_col=0)
polyadist = polyadist[[c for c in polyadist.columns if c not in EXCLUDE_SAMPLES]]

polyacounts = polyadist[polyadist.index >= MINIMUM_POLYA_LEN]

fig, ax = plt.subplots(1, 1, figsize=(7, 3.5))

xpositions = np.linspace(polyacounts.index[0],
                         polyacounts.index[-1], 100) ** XSCALE_TRANSFORM_FACTOR
xticks = xticks[(xticks >= polyacounts.index[0]) & (xticks <= polyacounts.index[-1])]

for (name, columns), color in zip(polyacounts.items(),
                                  colormap_lch(polyacounts.shape[1])):
    kde = gaussian_kde(columns.index ** XSCALE_TRANSFORM_FACTOR,
                       bw_method=KDE_BANDWIDTH, weights=columns.tolist())
    ax.plot(xpositions, kde(xpositions) * 100, c=color, label=name)

apply_dropped_spine(ax, xgrid=True)

ax.set_xlabel('Poly(A) length (nt)')
ax.set_ylabel('Density (a.u.)')
ax.set_xticks(xticks ** XSCALE_TRANSFORM_FACTOR)
ax.set_xticklabels(xticks)
ax.legend(fontsize=10, loc='best')

plt.subplots_adjust(bottom=.22, right=.95)

plt.savefig(snakemake.output[0])
