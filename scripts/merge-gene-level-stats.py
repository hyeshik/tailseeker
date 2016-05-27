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
import feather
import os

sm = snakemake

GENELEVEL_STATS_COLUMN_ORDER = {name: format(i, '08d') for i, name in enumerate("""
    polyA_mean polyA_mean_ci_lo polyA_mean_ci_hi polyA_median polyA_tag_count
""".split())}

GENEINFO_COLUMNS = "gene_name gene_biotype gene_description".split()

def load_sample_data(sample, filename):
    tbl = pd.read_csv(filename, index_col=0)
    tbl.columns = ['{}:{}'.format(sample, c) for c in tbl.columns]
    return tbl

def column_sortkey(name):
    sample, c = name.split(':', 1)
    return '{}:{}'.format(GENELEVEL_STATS_COLUMN_ORDER.get(c, c), sample)

# Load all tables
samples = sm.params.samples_by_genome[sm.wildcards.genome]
alldatatbl = pd.concat(list(map(load_sample_data, samples, sm.input)), axis=1)

# Merge selected columns from gene annotations
genes = feather.read_dataframe(os.path.join(sm.params.genomedir, 'annotations-gene.feather'))
genes = genes.set_index('gene_id')[GENEINFO_COLUMNS]

# Determine orders of columns and rows
cols_order = sorted(alldatatbl.columns, key=column_sortkey)
alldatatbl['__total_tags'] = alldatatbl[[c for c in alldatatbl.columns
                                        if c.endswith(':polyA_tag_count')]].sum(axis=1)
sorteddatatbl = alldatatbl.sort_values(by='__total_tags', ascending=False)[cols_order]
finaltbl = pd.concat([genes, sorteddatatbl], join='inner', axis=1)
del sorteddatatbl, alldatatbl, genes

# Write out in the csv format
finaltbl.to_csv(sm.output.csv)

# Write out in the feather format
feather.write_dataframe(finaltbl, sm.output.feather)
