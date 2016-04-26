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
import feather

stable_id_event_columns = """
old_stable_id old_version new_stable_id new_version
mapping_session_id type score""".split()

sie = pd.read_table(snakemake.input[0], compression='gzip',
                    names=stable_id_event_columns, na_values=['\\N'], keep_default_na=False,
                    low_memory=False)

sie_changes = sie.dropna(subset=['old_stable_id', 'new_stable_id']).copy()
sie_changes['old_id'] = sie_changes.apply(lambda r:
                            '{}.{}'.format(r['old_stable_id'], r['old_version']), axis=1)
sie_changes['new_id'] = sie_changes.apply(lambda r:
                            '{}.{}'.format(r['new_stable_id'], r['new_version']), axis=1)

# Prepare the initial round of old-new pairs discovery.
final_ids = set(sie_changes['new_id']) - set(sie_changes['old_id'])
is_final_mapping = sie_changes['new_id'].isin(final_ids)
final_mappings = [sie_changes[is_final_mapping][['old_id', 'new_id']]]
mapped = set(final_mappings[0]['new_id']) | set(final_mappings[0]['new_id'])

# Recursively resolve the eventual pair of old and new IDs.
for i in range(100):
    nested_mappings = pd.merge(final_mappings[-1],
                               sie_changes, left_on='old_id', right_on='new_id',
                               suffixes=['_new', '_old'])[['old_id_old', 'new_id_new']]
    nested_mappings.columns = ['old_id', 'new_id']
    nested_mappings = nested_mappings.drop_duplicates()
    nested_mappings = nested_mappings[~nested_mappings['old_id'].isin(mapped)]
    if len(nested_mappings) <= 0:
        break
    mapped |= set(nested_mappings['old_id'])
    final_mappings.append(nested_mappings)

# Finalize the mapping table.
fm = pd.concat(final_mappings, axis=0)
fm_with_type = (pd.merge(fm, sie_changes[['new_id', 'type']])
                  .drop_duplicates()
                  .sort_values(by=['old_id', 'new_id'])
                  .reset_index(drop=True))
feather.write_dataframe(fm_with_type, snakemake.output.versioned_table)

fm_with_type['old_id_simple'] = fm_with_type['old_id'].apply(lambda x: x.split('.', 1)[0])
fm_with_type['new_id_simple'] = fm_with_type['new_id'].apply(lambda x: x.split('.', 1)[0])
fm_simpler = fm_with_type[fm_with_type['old_id_simple'] != fm_with_type['new_id_simple']]
fm_simpler = fm_simpler[['old_id_simple', 'new_id_simple', 'type']].drop_duplicates()
fm_simpler.columns = ['old_id', 'new_id', 'type']
fm_simpler = fm_simpler.sort_values(by=['old_id', 'new_id']).reset_index(drop=True)
feather.write_dataframe(fm_simpler, snakemake.output.nonversioned_table)
