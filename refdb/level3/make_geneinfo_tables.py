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

from functools import reduce
import pandas as pd
import numpy as np
import feather
import re


def break_attribute_tags(attrs):
    alltags = reduce(lambda x, y: x | y, attrs.apply(lambda x: set(k for k, v in x)))
    cols = {}

    for tagname in alltags:
        # check if multiple tags can be included in a single annotation
        tags = attrs.apply(lambda x: [v for k, v in x if k == tagname])
        num_tags = tags.apply(len)
        dupcounts = num_tags.value_counts()
        is_plural = np.any(dupcounts.index > 1)
        is_required = (0 not in dupcounts.index)
        is_numeric = all(tags.apply(lambda x: not x or all(v.isdigit() for v in x)))

        if not is_numeric:
            if not is_plural:
                tags = tags.apply(lambda x: x[0] if x else np.nan)
            else:
                tags = tags.apply(lambda x: '; '.join(x))
        elif is_plural:
            tags = tags.apply(lambda x: list(map(int, x)))
        elif is_required:
            tags = tags.apply(lambda x: np.int64(x[0]))
            valuemax = tags.abs().max()
            if valuemax < 2**15:
                tags = tags.astype(np.int16)
            elif valuemax < 2**31:
                tags = tags.astype(np.int32)
            else:
                tags = tags.astype(np.int64)
        else:
            tags = tags.apply(lambda x: np.float64(x[0]) if x else np.nan)
            valuemax = tags.abs().max()
            if valuemax < 100000:
                tags = tags.astype(np.float32)
            else:
                tags = tags.astype(np.float64)

        cols[tagname] = tags

    return pd.DataFrame(cols)


def adopt_gtf_table(tbl):
    re_attrtags = re.compile('([^ ]+) "([^"]+)";')

    subtables = {}
    for feat_type, tagset in tbl.items():
        print(" - Inferring types for \"{}\"...".format(feat_type))
        attrs = tagset['attribute'].apply(lambda x: re_attrtags.findall(x))
        suppl = break_attribute_tags(attrs)
        subtables[feat_type] = pd.concat([tagset, suppl], axis=1)
        del subtables[feat_type]['attribute']

    return subtables



# == Load tables

print("Loading ENSEMBL tables...")
ensembl_tables = {}
tblreader = pd.read_table(snakemake.input.gtf, comment='#', compression='gzip', sep='\t',
                    names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
                           'fname', 'attribute'],
                    dtype={'seqname': str, 'source': str, 'feature': str,
                           'start': np.int32, 'end': np.int32, 'score': np.float32,
                           'strand': str,
                           'frame': np.float32, 'attribute': str},
                    na_values=['.'], keep_default_na=False, chunksize=65536,
                    low_memory=False)

for chunk in tblreader:
    for feat_type in chunk['feature'].unique():
        ensembl_tables.setdefault(feat_type, [])
        ensembl_tables[feat_type].append(chunk[chunk['feature'] == feat_type].copy())

for feat_type in list(ensembl_tables.keys()):
    ensembl_tables[feat_type] = pd.concat(ensembl_tables[feat_type], axis=0)


genedb_tbl_columns = """
    gene_id biotype analysis_id seq_region_id seq_region_start seq_region_end
    seq_region_strand display_xref_id source status description is_current
    canonical_transcript_id stable_id version created_date modified_date""".split()
genedb_tbl = pd.read_table(snakemake.input.database, compression='gzip',
                           names=genedb_tbl_columns, low_memory=False)


# == Make the gene info table from the GTF file.
print("Making a table from the GTF file...")
rtbl = adopt_gtf_table(ensembl_tables)

# == Split the UTR table into 5' and 3' UTRs if needed.
if 'UTR' in rtbl:
    cds_light = rtbl['CDS'][['strand', 'start', 'end', 'feature', 'transcript_id']]
    utr_light = rtbl['UTR'][['strand', 'start', 'end', 'feature', 'transcript_id']]
    cdsandutr = pd.concat([cds_light, utr_light], axis=0)

    utr5ix, utr3ix = [], []
    for trid, feats in cdsandutr.groupby('transcript_id'):
        cdsfeats = feats[feats['feature'] == 'CDS']
        utrfeats = feats[feats['feature'] == 'UTR']

        cds_start = cdsfeats['start'].min()
        cds_end = cdsfeats['end'].min()

        leftutrs = utrfeats[utrfeats['end'] < cds_start]
        rightutrs = utrfeats[utrfeats['start'] > cds_end]
        assert len(leftutrs.index & rightutrs.index) == 0
        assert len(feats['strand'].unique()) == 1

        if feats['strand'].iloc[0] == '+':
            utrleft, utrright = utr5ix, utr3ix
        else:
            utrleft, utrright = utr3ix, utr5ix

        utrleft.extend(leftutrs.index)
        utrright.extend(rightutrs.index)

    rtbl['five_prime_utr'] = cdsandutr.ix[utr5ix]
    rtbl['three_prime_utr'] = cdsandutr.ix[utr3ix]
    rtbl['five_prime_utr']['feature'] = 'five_prime_utr'
    rtbl['three_prime_utr']['feature'] = 'three_prime_utr'
    del rtbl['UTR']

for projtarget in ['gene', 'transcript']:
    # Check if gene_id in rtbl has version suffix.
    has_version_suffix = '.' in rtbl[projtarget]['gene_id'].iloc[0]
    if has_version_suffix:
        genedb_tbl['gene_id'] = genedb_tbl.apply(lambda row:
                                    '{}.{}'.format(row['stable_id'], row['version']), axis=1)
    else:
        genedb_tbl['gene_id'] = genedb_tbl['stable_id']

    rtbl[projtarget] = pd.merge(rtbl[projtarget], genedb_tbl[['gene_id', 'description']],
                                left_on='gene_id', right_on='gene_id', how='left')
    rtbl[projtarget]['gene_description'] = rtbl[projtarget]['description']
    del rtbl[projtarget]['description']


print("Saving tables to the disk...")
keycandidates = 'gene_id transcript_id seqname start end'.split()
for annotype, tbl in rtbl.items():
    if hasattr(snakemake.output, annotype):
        try:
            keys = [c for c in keycandidates if c in tbl.columns]
            tbl = tbl.sort_values(by=keys).reset_index(drop=True)
            feather.write_dataframe(tbl, getattr(snakemake.output, annotype))
        except:
            print(" - Writing feather failed! - {}".format(annotype))
            tbl.to_pickle(getattr(snakemake.output, annotype) + '.pkl')
            raise ValueError("Could not dump the gene info table using feather.")
    else:
        print(" - Skipping annotation type \"{}\".".format(annotype))

