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

import pandas as pd
import numpy as np
import lzma
from scipy.stats import t as stats_t
from tailseeker.stats import weighted_median
from operator import itemgetter


class TailGroupSummarizer:

    POLYA_WINDOWS = [ # 0-based, right non-inclusive
        [0, 1],
        [1, 4],
        [1, 6],
        [5, 11],
        [5, 16],
        [5, 26],
        [8, 16],
        [8, 26],
        [16, 26],
        [26, 41],
        [41, 100],
    ]

    def __init__(self, ci_width, polya_minimum_length=5):
        self.ci_width = ci_width
        self.polya_min_len = polya_minimum_length
        self.t = {}

    def mean_pa_length(self, polya_cnts):
        pa_cnts = polya_cnts[(polya_cnts.index >= self.polya_min_len) & (polya_cnts > 0)]

        n = pa_cnts.sum()
        if n == 0:
            return (np.nan, n, np.nan, np.nan)
        elif n == 1:
            return (pa_cnts.index[0], n, np.nan, np.nan)

        loglen = np.log(pa_cnts.index)
        pa_mean = np.average(loglen, weights=pa_cnts)
        pa_sem = (np.average((loglen - pa_mean) ** 2, weights=pa_cnts) / n) ** .5

        if (n - 1) not in self.t:
            self.t[n - 1] = stats_t.ppf(self.ci_width, df=n - 1)
        t = self.t[n - 1]

        pa_error = t * pa_sem
        return np.exp(pa_mean), n, np.exp(pa_mean - pa_error), np.exp(pa_mean + pa_error)

    def polya_length_stats(self, tails):
        cnt_by_pa = tails.sum(axis=2)
        mean_pa_packed = cnt_by_pa.apply(self.mean_pa_length)
        return pd.DataFrame({
            'polyA_mean': mean_pa_packed.apply(itemgetter(0)),
            'polyA_tag_count': mean_pa_packed.apply(itemgetter(1)),
            'polyA_mean_ci_lo': mean_pa_packed.apply(itemgetter(2)),
            'polyA_mean_ci_hi': mean_pa_packed.apply(itemgetter(3)),
            'polyA_median': cnt_by_pa.apply(weighted_median)
        })

    def polya_count_stats(self, tails):
        coldata = {}
        for left, right in self.POLYA_WINDOWS:
            modcounts = tails.ix[:, range(left, right)].sum(axis=1)
            coldata['polyA_tag_count_{}-{}'.format(left, right-1)] = modcounts.sum(axis=0)
        return pd.DataFrame(coldata)

    def polya_mods_stats(self, tails):
        coldata = {}
        for left, right in self.POLYA_WINDOWS:
            modcounts = tails.ix[:, range(left, right)].sum(axis=1)
            modsum = modcounts.multiply(np.array(modcounts.index), axis=0).sum(axis=0)
            avgmod = modsum / modcounts.sum(axis=0)
            coldata['average_mods_{}-{}'.format(left, right-1)] = avgmod

        return pd.DataFrame(coldata)

    def all_stats(self, tails):
        pastats = self.polya_length_stats(tails)
        modstats = self.polya_mods_stats(tails)
        return pd.concat([pastats, modstats], axis=1)



# ------------

confidence_interval = snakemake.params.confidence_interval_span

tgsum = TailGroupSummarizer(confidence_interval)
tailcounts = {
    'U': pd.read_msgpack(lzma.open(snakemake.input.U, 'rb')),
    'G': pd.read_msgpack(lzma.open(snakemake.input.G, 'rb')),
    'C': pd.read_msgpack(lzma.open(snakemake.input.C, 'rb')),
}

stats = [
    tgsum.polya_length_stats(tailcounts['U']),
    tgsum.polya_count_stats(tailcounts['U']),
]

for modtype, tailcnt in tailcounts.items():
    tbl = tgsum.polya_mods_stats(tailcnt)
    tbl.columns = [col.replace('_mods_', '_{}_'.format(col)) for col in tbl.columns]
    stats.append(tbl)

(pd.concat(stats, axis=1)
   .sort_values(by='polyA_tag_count', ascending=False)
   .to_csv(snakemake.output[0]))

