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

from tailseeker import fileutils
import subprocess as sp
import pandas as pd
import numpy as np
import os.path
from io import TextIOWrapper
from textwrap import wrap
import feather

SEQLINEWIDTH = 72
SEQFLANKLENGTH = 200

# Calculate lengths of 5'UTRs and 3'UTRs
fivep_utr = feather.read_dataframe(snakemake.input.fivep_utr)
fivep_utr['length'] = fivep_utr['end'] - fivep_utr['start'] + 1
fivep_utr_length = fivep_utr.groupby('transcript_id').agg({'length': 'sum'})['length']
threep_utr = feather.read_dataframe(snakemake.input.threep_utr)
threep_utr['length'] = threep_utr['end'] - threep_utr['start'] + 1
threep_utr_length = threep_utr.groupby('transcript_id').agg({'length': 'sum'})['length']
utrlengths = pd.DataFrame({'fivep_utr_length': fivep_utr_length,
                           'threep_utr_length': threep_utr_length}).fillna(0)
del fivep_utr, threep_utr, fivep_utr_length, threep_utr_length

genes = feather.read_dataframe(snakemake.input.gene)
exons = feather.read_dataframe(snakemake.input.exon)

with fileutils.TemporaryDirectory() as tmpdir:
    tmpbed = os.path.join(tmpdir, 'exons.bed')

    print('\n'.join(exons.sort_values(by=['seqname', 'start'])
                    .apply(lambda row: '{}\t{}\t{}\t{}-{}\t.\t{}'.format(
                            row['seqname'], row['start'] - 1, row['end'],
                            row['transcript_id'], row['exon_number'],
                            row['strand']), axis=1)),
          file=open(tmpbed, 'w'))

    tmpfasta = os.path.join(tmpdir, 'exons.fa')
    sp.check_call(['bedtools', 'getfasta', '-fi', snakemake.input.genome,
                   '-bed', tmpbed, '-tab', '-s', '-name', '-fo',
                   tmpfasta])

    exonseqs = pd.read_table(tmpfasta, names=['exon_name', 'sequence'])

exonseqs['transcript_id'] = exonseqs['exon_name'].apply(lambda x: x.split('-')[0])
exonseqs['exon_number'] = exonseqs['exon_name'].apply(lambda x: int(x.split('-')[1]))
exonseqs = exonseqs.sort_values(by=['transcript_id', 'exon_number'])

def alter_UTR_cases(row):
    fivep_utr = (row['sequence'][:int(row['fivep_utr_length'])]
                 if row['fivep_utr_length'] > 0 else '').lower()
    threep_utr = (row['sequence'][-int(row['threep_utr_length']):]
                 if row['threep_utr_length'] > 0 else '').lower()
    cds = row['sequence'][len(fivep_utr):len(row['sequence'])-len(threep_utr)]
    fullseq = fivep_utr + cds + threep_utr
    assert len(fullseq) == len(row['sequence'])
    return fullseq

transcriptseqs = (
    exonseqs.groupby('transcript_id')
            .agg({'sequence': 'sum'})
            .join(utrlengths)
            .reset_index())
transcriptseqs['seqaugmented'] = transcriptseqs.apply(alter_UTR_cases, axis=1)

gzwriter = sp.Popen(['bgzip', '-c'], stdin=sp.PIPE, stdout=open(snakemake.output.fasta, 'wb'))
print(''.join(transcriptseqs.apply(lambda row:
        '>{}\n{}\n'.format(row['transcript_id'],
                           '\n'.join(wrap(row['seqaugmented'], SEQLINEWIDTH))), axis=1)),
      file=TextIOWrapper(gzwriter.stdin))

print(''.join(transcriptseqs.apply(lambda row:
        '{}\t{}\n'.format(row['transcript_id'], row['seqaugmented'][-SEQFLANKLENGTH:]), axis=1)),
      file=open(snakemake.output.threependtable, 'w'))

