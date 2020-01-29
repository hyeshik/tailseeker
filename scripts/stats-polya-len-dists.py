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

from tailseeker.tabledefs import taginfo, refined_taginfo
from io import BytesIO
import subprocess as sp
import pandas as pd
import numpy as np


def get_polya_length_hist(pacallsfile, badflagmask, maxpalen, columnno):
    freqoutput = sp.check_output(
        'zcat {input_file} | perl -e \'while (<>) {{ \
            chomp; my @f=split; if ((@f[2] & {flagmask}) == 0) {{ \
                print "$f[{polya_column}]\\n"; }} }}\' | \
            sort -n | uniq -c'.format(
            input_file=pacallsfile, polya_column=columnno,
            flagmask=badflagmask), shell=True)

    pacounts = pd.read_table(BytesIO(freqoutput), skipinitialspace=True, sep=' ',
                             names=['count', 'polya'], dtype=np.int64).set_index('polya')

    pacounts = pacounts.reindex(list(range(maxpalen + 1)), fill_value=0)

    return pacounts['count']


# ===

inputfiles = list(zip(snakemake.params.samplenames, snakemake.input))

tableschema = refined_taginfo if snakemake.params.refined else taginfo
columnname = 'unaligned_polyA' if snakemake.params.refined else 'polyA'
columnno = tableschema['names'].index(columnname)

pacounts = pd.DataFrame()
for samplename, lencallfile in inputfiles:
    pacounts[samplename] = get_polya_length_hist(
        lencallfile, snakemake.params.badflagmask, snakemake.params.maxpalength, columnno)

# trim all-zero rows at the end of the list
effective_maximum_pa_len = np.where(pacounts.sum(axis=1) > 0)[0].max()
pacounts = pacounts[pacounts.index <= effective_maximum_pa_len]

# write the output in csv format
pacounts.to_csv(snakemake.output[0])

