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

TARGETS.extend(expand('associations/{sample}.txt.gz', sample=EXP_SAMPLES))

rule associate_tags_to_genes:
    input: 'alignments/{sample}_single.bam'
    output: temp('scratch/associations/{sample}.txt')
    threads: 4
    run:
        genomedir = os.path.join(TAILSEEKER_DIR, 'refdb', 'level3',
                                 CONF['reference_set'][wildcards.sample])
        required_mapq = CONF['gene_level_stats']['required_mapping_quality']

        shell('{SAMTOOLS_CMD} view -b -F 260 -q {required_mapq} {input} | \
               {BEDTOOLS_CMD} intersect -nonamecheck -split -wa -wb -bed \
                -a - -b {genomedir}/exons.gtf.gz -s | \
               cut -f4,21 | \
               {PARALLEL_CMD} -j {threads} --pipe --block 10M \
                    "sed -e \'s,gene_id \\"\\([^\\"]*\\)\\".*$,\\1,g\'" | \
               sort -k1,2 | uniq > {output}')

rule count_multiple_associations:
    input: 'scratch/associations/{sample}.txt'
    output: 'associations/{sample}.txt.gz'
    threads: 3
    shell: 'cut -d\'\t\' -f1 {input} | uniq -c | \
            awk \'BEGIN {{ OFS="\t"; }} {{ print $2, $1; }}\' | \
            join -j 1 -t\'\t\' {input} /dev/stdin | \
            {PARALLEL_CMD} -j {threads} --pipe --block 10M \
                "sed -e \'s,^\\([^:]*\\):0*\\([0-9][0-9]*\\):[^\t]*\t,\\1\t\\2\t,g\'" | \
            {BGZIP_CMD} -@ {threads} -c > {output}'


TARGETS.extend(expand('tagcounts/{sample}-{ambigtype}-{modtype}.msgpack.xz',
                      sample=EXP_SAMPLES, ambigtype=['single', 'multi'],
                      modtype=['U', 'C', 'G']))

rule make_gene_level_counts:
    input:
        taginfo='refined-taginfo/{sample}.txt.gz',
        associations='associations/{sample}.txt.gz'
    output: 'tagcounts/{sample}-{ambigtype}-{modtype}.msgpack.xz'
    run:
        from tailseeker import tabledefs
        import pandas as pd
        import numpy as np
        import lzma

        assoctbl = pd.read_table(input.associations, **tabledefs.associations)
        taginfotbl = pd.read_table(input.taginfo, **tabledefs.refined_taginfo)

        tbl = pd.merge(assoctbl, taginfotbl, how='inner', left_on=['tile', 'cluster'],
                       right_on=['tile', 'cluster'])

        bad_flags = CONF['gene_level_stats']['bad_flags_filter']
        max_modcount = CONF['gene_level_stats']['maximum_nonA_mod_count']
        delim_settings = CONF['delimiter'][wildcards.sample]
        min_preamble_length = delim_settings[0] - 1 + len(delim_settings[1]) - 1
        R3 = CONF['read_cycles']['R3']
        max_polya = R3[1] - R3[0] + 1 - min_preamble_length
        mod_of_interest = wildcards.modtype

        tbl[mod_of_interest] = tbl[mod_of_interest].clip_upper(max_modcount)
        if wildcards.ambigtype == 'single':
            tbl = tbl[tbl['ambig'] <= 1]
        filtered_tbl = tbl[(tbl['pflags'] & bad_flags) == 0]

        polya_axis = np.arange(max_polya + 1).astype(np.int32)
        modcount_axis = np.arange(max_modcount + 1).astype(np.int32)

        counts_dfs = {}
        for geneid, tags in filtered_tbl.groupby('gene'):
            tagcounts = tags.groupby(['polyA', mod_of_interest]).agg('count').reset_index()
            countsgrid = (pd.pivot_table(tagcounts, index='polyA', columns=mod_of_interest,
                                          values='cluster', fill_value=0).astype(np.uint32)
                                   .reindex_axis(polya_axis, axis=0, fill_value=0)
                                   .reindex_axis(modcount_axis, axis=1, fill_value=0))
            counts_dfs[geneid] = countsgrid

        counts_dfs = pd.Panel(counts_dfs)
        counts_dfs.to_msgpack(lzma.open(output[0], 'wb'))

# ex: syntax=snakemake
