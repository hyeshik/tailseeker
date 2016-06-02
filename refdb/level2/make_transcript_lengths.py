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

import pandas as pd
import numpy as np
import re

gtf_input = snakemake.input[0]

gtf = pd.read_table(gtf_input, comment='#', compression='gzip', sep='\t',
                    names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
                           'fname', 'attribute'],
                    dtype={'seqname': str, 'source': str, 'feature': str,
                           'start': np.int32, 'end': np.int32, 'score': np.float32,
                           'strand': str,
                           'frame': np.float32, 'attribute': str},
                    na_values=['.'], keep_default_na=False, chunksize=65536,
                    low_memory=False)

pat_transcript_id = re.compile('transcript_id "([^"]+)"')

exons = []
for chunk in gtf:
    chunkexons = chunk[chunk['feature'] == 'exon'].copy()
    chunkexons['transcript_id'] = chunkexons['attribute'].apply(lambda x:
                                                        pat_transcript_id.findall(x)[0])
    del chunkexons['attribute']
    exons.append(chunkexons)

exons = pd.concat(exons, axis=0)
exons['length'] = exons['end'] - exons['start'] + 1
exons.groupby('transcript_id').agg({'length': 'sum'}).to_csv(
        snakemake.output[0], sep='\t', header=False)

