#!/usr/bin/env python3
#
# Copyright (c) 2016 Institute for Basic Science
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

from tailor.powersnake import *
import pandas as pd


def load_tile_demux_table(tile, filename):
    tbl = pd.read_csv(filename, na_values=None, keep_default_na=False)
    tbl.insert(0, 'tile', tile)
    return tbl


def process(tile_table_pairs, output):
    tilestats = pd.concat(
        [load_tile_demux_table(tile, tablefile)
         for tile, tablefile in tile_table_pairs], axis=0)

    SUM_FIELDS = 'cln_no_index_mm cln_1_index_mm cln_2+_index_mm cln_no_delim'.split()
    aggfuncs = {c: 'sum' for c in SUM_FIELDS}
    aggfuncs.update({c: 'first' for c in tilestats.columns
                     if c not in SUM_FIELDS and c != 'name'})

    totalsum = tilestats.groupby('name').agg(aggfuncs).reset_index()[tilestats.columns]
    totalsum['tile'] = '(total)'

    merged = pd.concat([totalsum,
                        tilestats.sort_values(by=['name', 'tile'])], axis=0)
    merged['cln_passed'] = (
        merged['cln_no_index_mm cln_1_index_mm cln_2+_index_mm'.split()].sum(axis=1) -
        merged['cln_no_delim'])
    merged.to_csv(output, index=False)


## XXX Implement the standalone interface.

if is_snakemake_child:
    inputfiles = list(zip(params.tiles, input))
    process(inputfiles, output[0])
