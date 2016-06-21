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

GENEINFO_COLUMNS = "gene_name gene_type gene_description".split()
GENETYPE_ALIASES = ['gene_biotype']

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
for alias in GENETYPE_ALIASES:
    if alias in genes.columns:
        genes['gene_type'] = genes[alias]
genes = genes.set_index('gene_id')[GENEINFO_COLUMNS]

# Determine orders of columns and rows
cols_order = GENEINFO_COLUMNS + sorted(alldatatbl.columns, key=column_sortkey)
alldatatbl['__total_tags'] = alldatatbl[[c for c in alldatatbl.columns
                                        if c.endswith(':polyA_tag_count')]].sum(axis=1)
finaltbl = pd.merge(genes, alldatatbl, how='right', left_index=True,
                    right_index=True).sort_values(
                    by='__total_tags', ascending=False)[cols_order]
del alldatatbl, genes

# Write out in the csv format
print("Writing to " + sm.output.csv)
finaltbl.to_csv(sm.output.csv)

# Write out in the feather format
print("Writing to " + sm.output.feather)
feather.write_dataframe(finaltbl.reset_index(), sm.output.feather)

# Write out in the Excel format
EXCEL_FORMAT_OPTS = {
    'gene_name': [10.5, {'bg_color': '#ededed'}],
    'gene_type': [8.5, {'font_color': '#787878'}],
    'gene_description': [30, {'align': 'left'}],
    ':polyA_mean': [6.2,
        {'align': 'right', 'bg_color': '#ddebf7', 'num_format': '0.0'}],
    ':polyA_mean_ci_lo': [6.2,
        {'align': 'right', 'bg_color': '#fff2cc', 'num_format': '0.0'}],
    ':polyA_mean_ci_hi': [6.2,
        {'align': 'right', 'bg_color': '#fce4d6', 'num_format': '0.0'}],
    ':polyA_median': [6.2,
        {'align': 'right', 'num_format': '0.0'}],
    ':polyA_tag_count': [6.2,
        {'align': 'right', 'bg_color': '#e2efda', 'num_format': '0.0'}],
}
EXCEL_HEADER_FORMAT = {'align': 'left', 'text_wrap': True, 'bold': True,
                       'valign': 'top', 'bg_color': '#d9d9d9'}
EXCEL_INDEX_FORMAT = {'align': 'left', 'bold': False}

with pd.ExcelWriter(sm.output.excel, engine='xlsxwriter') as writer:
    print("Preparing an Excel table.")
    finaltbl.to_excel(writer, sheet_name='TAIL-seq')
    book = writer.book
    sheet = writer.sheets['TAIL-seq']

    print("Prettify'ing the table.")
    for colname, (width, formatting) in EXCEL_FORMAT_OPTS.items():
        if colname == '':
            colnums = [0]
        else:
            colnums = [coli+1 for coli, col in enumerate(finaltbl.columns)
                       if col.endswith(colname)]

        fmt = book.add_format(formatting)
        sheet.set_column(colnums[0], colnums[-1], width, fmt)

    for i, name in enumerate(['gene_id'] + list(finaltbl.columns)):
        sheet.write_string(0, i, name, book.add_format(EXCEL_HEADER_FORMAT))
    for i, name in enumerate(finaltbl.index, 1):
        sheet.write_string(i, 0, name, book.add_format(EXCEL_INDEX_FORMAT))
    sheet.set_column(0, 0, 18.5)

    sheet.freeze_panes('C2')
    sheet.set_zoom(115)

    print("Writing to " + sm.output.excel)
    writer.save()
