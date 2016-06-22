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

from tailseeker import tabledefs
import pandas as pd
import numpy as np
import lzma

sm = snakemake

print("Loading required tables...")
assoctbl = pd.read_table(sm.input.associations, **tabledefs.associations)
taginfotbl = pd.read_table(sm.input.taginfo, **tabledefs.refined_taginfo)

print("Joining tables...")
tbl = pd.merge(assoctbl, taginfotbl, how='inner', left_on=['tile', 'cluster'],
               right_on=['tile', 'cluster'])

min_preamble_length = sm.params.delim_settings[0] - 1 + len(sm.params.delim_settings[1]) - 1
max_polya = sm.params.R3[1] - sm.params.R3[0] + 1 - min_preamble_length
mod_of_interest = sm.wildcards.modtype

print("Filtering out bad tags...")
tbl[mod_of_interest] = tbl[mod_of_interest].clip_upper(sm.params.max_modcount)
if sm.wildcards.ambigtype == 'single':
    tbl = tbl[tbl['ambig'] <= 1]
filtered_tbl = tbl[(tbl['pflags'] & sm.params.bad_flags) == 0]
del tbl

def write_tagcounts(output, tbl, polya_axis, modcount_axis):
    counts_dfs = {}
    groups = tbl.groupby('gene')
    for i, (geneid, tags) in enumerate(groups):
        if i % 100 == 0:
            print(" - Processing genes... ({0}/{1})     \r".format(i, len(groups)), end='')
        tagcounts = tags.groupby(['polyA', mod_of_interest]).agg('count').reset_index()
        countsgrid = (pd.pivot_table(tagcounts, index='polyA', columns=mod_of_interest,
                                      values='cluster', fill_value=0).astype(np.uint32)
                               .reindex_axis(polya_axis, axis=0, fill_value=0)
                               .reindex_axis(modcount_axis, axis=1, fill_value=0))
        counts_dfs[geneid] = countsgrid

    print("\n - Writing to a file.")
    counts_dfs = pd.Panel(counts_dfs)
    counts_dfs.to_msgpack(lzma.open(output, 'wb'))


# Intact poly(A) tails

polya_axis = np.arange(max_polya + 1).astype(np.int32)
modcount_axis = np.arange(sm.params.max_modcount + 1).astype(np.int32)

print("Writing out the canonical tails table...")
is_intact_tail = ((filtered_tbl['polyA'] >= sm.params.polyA_assume_intact) |
                  ((filtered_tbl['pflags'] & tabledefs.PAFLAG_LIKELY_HAVE_INTACT_END) != 0))
write_tagcounts(sm.output.canonical, filtered_tbl[is_intact_tail], polya_axis, modcount_axis)


# Non-poly(A) or degraded tails

polya_axis = np.arange(sm.params.polyA_assume_intact).astype(np.int32)
modcount_axis = np.arange(sm.params.max_modcount + 1).astype(np.int32)

print("Writing out the non-canonical tails table...")
write_tagcounts(sm.output.noncanonical, filtered_tbl[~is_intact_tail], polya_axis, modcount_axis)
