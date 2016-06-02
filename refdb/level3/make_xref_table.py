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
import feather


column_names = {
    'gene': '''
        gene_id biotype analysis_id seq_region_id seq_region_start seq_region_end
        seq_region_strand display_xref_id source status description is_current
        canonical_transcript_id stable_id version created_date modified_date'''.split(),
    'transcript': '''
        transcript_id gene_id analysis_id seq_region_id seq_region_start seq_region_end
        seq_region_strand display_xref_id source biotype status description is_current
        canonical_translation_id stable_id version created_date modified_date'''.split(),
    'xref': '''
        xref_id external_db_id dbprimary_acc display_label version description
        info_type info_text'''.split(),
    'external_db': '''
        external_db_id db_name db_release status priority db_display_name type
        secondary_db_name secondary_db_table description'''.split(),
    'translation': '''
        translation_id transcript_id seq_start start_exon_id seq_end end_exon_id
        stable_id version created_date modified_date'''.split(),
    'identity_xref': '''
        object_xref_id xref_identity ensembl_identity xref_start xref_end
        ensembl_start ensembl_end cigar_line score evalue'''.split(),
    'object_xref': '''
        object_xref_id ensembl_id ensembl_object_type xref_id linkage_annotation
        analysis_id'''.split(),
    'ontology_xref': '''
        object_xref_id source_xref_id linkage_type'''.split(),
}


# == Load tables

print("Loading ENSEMBL tables...")
annotbls = {
    tblname: pd.read_table(snakemake.input[tblname], compression='gzip', low_memory=False,
                           names=columns, na_values=['\\N'], keep_default_na=False,
                           dtype={col: object for col in columns if col.endswith('_id')})
    for tblname, columns in column_names.items()}


# == Make the DNA-RNA-Protein mapping table.

print("Making the internal cross-refernece table...")
tmp = pd.merge(annotbls['gene'][['stable_id', 'gene_id']],
               annotbls['transcript'][['stable_id', 'transcript_id', 'gene_id']],
               left_on='gene_id', right_on='gene_id', suffixes=('_gene', '_transcript'))
ensembl_id_internal_refs = pd.merge(tmp,
                annotbls['translation'][['stable_id', 'transcript_id', 'translation_id']],
                left_on='transcript_id', right_on='transcript_id')[
                    ['stable_id_gene', 'stable_id_transcript', 'stable_id']]
ensembl_id_internal_refs.columns = ['ensembl_id_gene', 'ensembl_id_transcript',
                                    'ensembl_id_translation']
ensembl_id_internal_refs = ensembl_id_internal_refs.sort_values(
                                by=['ensembl_id_gene', 'ensembl_id_transcript',
                                    'ensembl_id_translation']).reset_index(drop=True)
feather.write_dataframe(ensembl_id_internal_refs, snakemake.output.internal_xrefs)


# == Prepare fundamental tables for xref tables.
print("Preparing fundamental tables for external references...")
xrefs = pd.merge(annotbls['external_db'][['external_db_id', 'db_name']],
                 annotbls['xref'][['xref_id', 'external_db_id', 'dbprimary_acc']])
xrefs.columns = ['external_db_id', 'xref_source', 'xref_id', 'xref_external_id']
xrefs = xrefs[['xref_source', 'xref_id', 'xref_external_id']].copy()

tmp = pd.merge(annotbls['object_xref'], annotbls['identity_xref'],
               left_on='object_xref_id', right_on='object_xref_id')
identity_xrefs = pd.merge(tmp, xrefs, left_on='xref_id', right_on='xref_id')

tmp = pd.merge(annotbls['object_xref'], annotbls['ontology_xref'],
               left_on='object_xref_id', right_on='object_xref_id')
ontology_xrefs = pd.merge(tmp, xrefs, left_on='xref_id', right_on='xref_id')


# == Make identity xrefs.
print("Making the identity cross-references...")
for grptype, grp in identity_xrefs.groupby('ensembl_object_type'):
    if grptype == 'Gene':
        gene_identity_xrefs = pd.merge(annotbls['gene'], grp,
                                       left_on='gene_id', right_on='ensembl_id',
                                       suffixes=['', '_idxr'])
    elif grptype == 'Transcript':
        transcript_identity_xrefs = pd.merge(annotbls['transcript'], grp,
                                             left_on='transcript_id', right_on='ensembl_id',
                                             suffixes=['', '_idxr'])
    elif grptype == 'Translation':
        translation_identity_xrefs = pd.merge(annotbls['translation'], grp,
                                              left_on='translation_id', right_on='ensembl_id',
                                              suffixes=['', '_idxr'])
    else:
        raise ValueError("Unknown object type {}".format(grptype))

# Broadcast xrefs to related objects (ie. from protein to gene).
direct_xrefs = [
    ('gene',        gene_identity_xrefs),
    ('transcript',  transcript_identity_xrefs),
    ('translation', translation_identity_xrefs),
]

identity_xrefs = []
for name_from, xrefs_from in direct_xrefs:
    dtbl = xrefs_from[['stable_id', 'xref_source', 'xref_external_id']].copy()
    dtbl.columns = ['ensembl_id', 'xref_source', 'xref_external_id']
    dtbl['ensembl_basis'] = name_from

    for name_to, xrefs_to in direct_xrefs:
        if name_from == name_to:
            identity_xrefs.append(dtbl)
            continue

        expanded = pd.merge(dtbl, ensembl_id_internal_refs[['ensembl_id_' + name_from,
                                                            'ensembl_id_' + name_to]],
                            left_on='ensembl_id', right_on='ensembl_id_' + name_from)
        expanded = expanded[['ensembl_id_' + name_to, 'xref_source', 'xref_external_id',
                             'ensembl_basis']].copy()
        expanded.columns = dtbl.columns
        identity_xrefs.append(expanded)

identity_xrefs = pd.concat(identity_xrefs, axis=0)

print("Saving the table to disk...")
identity_xrefs = identity_xrefs.sort_values(by=['ensembl_id', 'xref_source',
                                                'xref_external_id']).reset_index(drop=True)
feather.write_dataframe(identity_xrefs, snakemake.output.identity_xrefs)
del identity_xrefs


# == Make ontology xrefs.
print("Making the ontology cross-references...")
for grptype, grp in ontology_xrefs.groupby('ensembl_object_type'):
    if grptype == 'Translation':
        translation_ontology_xrefs = pd.merge(annotbls['translation'], grp,
                                              left_on='translation_id', right_on='ensembl_id',
                                       suffixes=['', '_onxr'])
    else:
        raise ValueError("Unknown object type {}".format(grptype))

go_xrefs_direct = translation_ontology_xrefs[['stable_id', 'xref_source', 'xref_external_id',
                                              'linkage_type']].copy()
go_xrefs_direct.columns = ['ensembl_id', 'xref_source', 'xref_external_id', 'linkage_type']
go_xrefs_direct['ensembl_basis'] = 'translation'

# Broadcast xrefs to related objects (ie. from protein to gene).
go_xrefs_ext = pd.merge(go_xrefs_direct, ensembl_id_internal_refs,
                        left_on='ensembl_id', right_on='ensembl_id_translation')

go_xrefs_for_gene = go_xrefs_ext[['ensembl_id_gene', 'xref_source', 'xref_external_id',
                                  'linkage_type', 'ensembl_basis']].copy().drop_duplicates()
go_xrefs_for_gene.columns = ['ensembl_id', 'xref_source', 'xref_external_id', 'linkage_type',
                             'ensembl_basis']

go_xrefs_for_tr = go_xrefs_ext[['ensembl_id_transcript', 'xref_source', 'xref_external_id',
                                'linkage_type', 'ensembl_basis']].copy().drop_duplicates()
go_xrefs_for_tr.columns = ['ensembl_id', 'xref_source', 'xref_external_id', 'linkage_type',
                           'ensembl_basis']

ontology_xrefs = pd.concat([go_xrefs_direct, go_xrefs_for_tr, go_xrefs_for_gene], axis=0)

print("Saving the table to disk...")
ontology_xrefs = ontology_xrefs.sort_values(by=['ensembl_id', 'xref_source',
                                                'xref_external_id']).reset_index(drop=True)
feather.write_dataframe(ontology_xrefs, snakemake.output.ontology_xrefs)

