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
from csv import QUOTE_NONE
from functools import reduce


def quote_internal_trid(trid):
    # Internal transcript ID in Xenbase includes ':' character and
    # even doesn't keep identical prefix (PAC or PAC4GC) in different features.
    return 'IBSXL{:010d}'.format(int(trid.split(':', 1)[1]))

def process(inputfile, outputfile):
    gff3 = pd.read_table(inputfile, comment='#',
                         names=['seqname', 'source', 'feature', 'start', 'end',
                                'score', 'strand', 'frame', 'attribute'])
    gff3['attrs'] = gff3['attribute'].apply(lambda x: dict(tok.split('=', 1)
                                                           for tok in x.split(';')))

    # Split features into feature type groups.
    byfeat = {}
    for featname, grp in gff3.groupby('feature'):
        grp = grp.copy()
        attrnames = reduce(lambda x, y: x|y, grp['attrs'].apply(set))
        for attrname in attrnames:
            grp['_' + attrname] = grp['attrs'].apply(lambda x: x.get(attrname, np.nan))
        del grp['attrs'], grp['attribute']
        byfeat[featname] = grp

    # Process genes
    byfeat['gene']['gene_id'] = byfeat['gene']['_ID']
    byfeat['gene']['gene_name'] = byfeat['gene']['_Name']
    byfeat['gene'].head()

    converted = {'gene': byfeat['gene']}

    # Process mRNAs
    mRNA_merged = pd.merge(byfeat['mRNA'], byfeat['gene'], how='left',
                           left_on='_Parent', right_on='_ID', suffixes=['', '_gene_'])
    mRNA_merged['feature'] = 'transcript'
    mRNA_merged['transcript_name'] = mRNA_merged['_Name']
    mRNA_merged['transcript_biotype'] = 'protein_coding'
    mRNA_merged['transcript_id'] = list(map(quote_internal_trid, mRNA_merged['_ID']))
    mRNA_merged['transcript_is_longest'] = mRNA_merged['_longest']
    mRNA_merged = mRNA_merged[[col for col in mRNA_merged.columns if not col.endswith('_')]]
    converted['transcript'] = mRNA_merged

    # Create exons from CDS and UTRs.
    exons = pd.concat([byfeat['CDS'], byfeat['five_prime_UTR'], byfeat['three_prime_UTR']])

    processed_exons = []
    for transcriptid, regions in exons.groupby('_Parent'):
        regions = regions.sort_values(by=['start', 'end'])
        mgexons = []
        for _, row in regions.iterrows():
            if not mgexons or mgexons[-1].end + 1 != row.start:
                mgexons.append(row)
            else:
                mgexons[-1].end = row.end

        mgexons = pd.DataFrame.from_records(mgexons)
        mgexons['feature'] = 'exon'
        mgexons['_ID'] = ['{}.exon.{}'.format(quote_internal_trid(transcriptid), i+1)
                          for i in range(len(mgexons))]
        processed_exons.append(mgexons)

    processed_exons = pd.concat(processed_exons, axis=0).reset_index(drop=True)
    processed_exons['_Parent'] = list(map(quote_internal_trid, processed_exons['_Parent']))

    exons_tr = pd.merge(processed_exons, converted['transcript'], how='left',
                        left_on='_Parent', right_on='transcript_id', suffixes=['', '_tr_'])
    exons_tr = exons_tr[[col for col in exons_tr.columns if not col.endswith('_')]]
    exons_tr['exon_id'] = exons_tr['_ID']
    converted['exon'] = exons_tr

    # Process CDS and UTRs
    for regtype in ['CDS', 'five_prime_UTR', 'three_prime_UTR']:
        byfeat[regtype]['_Parent'] = list(map(quote_internal_trid, byfeat[regtype]['_Parent']))
        reg_tr = pd.merge(byfeat[regtype], converted['transcript'], how='left',
                          left_on='_Parent', right_on='transcript_id',
                          suffixes=['', '_tr_'])
        reg_tr = reg_tr[[col for col in reg_tr.columns if not col.endswith('_')]]
        regtype_ensembl = regtype if regtype == 'CDS' else regtype.lower()
        reg_tr['feature'] = regtype_ensembl
        converted[regtype_ensembl] = reg_tr

    # Make the ENSEMBL-like GTF rows.
    def encode_attributes(tbl):
        attrs = tbl[[col for col in tbl.columns[8:] if not col.startswith('_')]]
        encoded_attrs = attrs.apply(lambda row: ' '.join(
                            '{} "{}";'.format(name, value.replace(' ', ''))
                            for name, value in row.items()), axis=1)
        ret = tbl.iloc[:, :8].copy()
        ret['attribute'] = encoded_attrs
        return ret

    feats_encoded = pd.concat([encode_attributes(feats)
                               for feats in converted.values()], axis=0)
    feats_encoded = feats_encoded.sort_values(by=['seqname', 'start', 'end', 'feature'])
    feats_encoded.to_csv(outputfile, sep='\t', index=False, header=False,
                         doublequote=False, quotechar="'", quoting=QUOTE_NONE)


if __name__ == '__main__':
    import sys
    process(sys.stdin, sys.stdout)
