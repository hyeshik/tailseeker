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


def prepare_experiment_groups_view():
    groups = defaultdict(list)
    for sample, group in CONF['experiment_groups'].items():
        groups[str(group)].append(sample)
    return dict(groups)
EXPERIMENT_GROUPS = prepare_experiment_groups_view()


# ---
# Sequence alignment to contaminants (optional).
# ---

rule contaminant_alignment:
    input: 'fastq/{sample}_R5.fastq.gz'
    output: temp('scratch/contaminants-unmapped/{sample}.txt')
    threads: THREADS_MAXIMUM_CORE
    params: scratch='scratch/contaminants-{sample}'
    run:
        genomedir = os.path.join(TAILSEEKER_DIR, 'refdb', 'level2',
                        CONF['reference_set'][wildcards.sample], 'contaminants.star')

        if os.path.isdir(params.scratch):
            shutil.rmtree(params.scratch)
        os.makedirs(params.scratch)

        shell('{STAR_CMD} --runThreadN {threads} \
                --genomeDir {genomedir} --outSAMtype None --outSAMmode None \
                --outFilterMultimapNmax 512 --outReadsUnmapped Fastx \
                --readFilesIn {input} --readFilesCommand zcat \
                --outTmpDir {params.scratch}/tmp \
                --outFileNamePrefix {params.scratch}/')

        shell('split -n r/1/4 {params.scratch}/Unmapped.out.mate1 | \
               colrm 1 1 | sort > {output}')
        shutil.rmtree(params.scratch)


if CONF['read_filtering']['contaminant_filtering']:
    TARGETS.extend(expand('fastq-filtered/{sample}_{read}.fastq.gz',
                          sample=EXP_SAMPLES, read=['R5', 'R3']))

rule make_filtered_fastq:
    input:
        fastq='fastq/{sample}_{read}.fastq.gz',
        survivorids='scratch/contaminants-unmapped/{sample}.txt'
    output: 'fastq-filtered/{sample}_{read,[^_]+}.fastq.gz'
    threads: 2
    shell: '{SEQTK_CMD} subseq {input.fastq} {input.survivorids} | \
            {BGZIP_CMD} -@ {threads} -c > {output}'

# ---
# Sequence alignment to the genome.
# ---

TARGETS.extend(expand('alignments/{sample}_{type}.bam{suffix}',
                      sample=EXP_SAMPLES, type=['single', 'paired'], suffix=['', '.bai']))

def inputs_for_STAR_alignment(wildcards):
    fastq_filename = 'fastq{filtered_suffix}/{sample}_{{read}}.fastq.gz'.format(
                            sample=wildcards.sample,
                            filtered_suffix='-filtered'
                                            if CONF['read_filtering']['contaminant_filtering']
                                            else '')

    inputs = [fastq_filename.format(read='R5')]
    if wildcards.type == 'paired':
        inputs.append(fastq_filename.format(read='R3'))
    return inputs

rule STAR_alignment:
    input: inputs_for_STAR_alignment
    output:
        mapped=temp('scratch/alignments/{sample}_STAR_{type,[^_.]+}.bam'),
        unmapped5=temp('scratch/unmapped-reads/{sample}-{type}_R5.fastq.gz'),
        unmapped3=temp('scratch/unmapped-reads/{sample}-{type}_R3.fastq.gz'),
        transcriptome=temp('scratch/tr-alignments/{sample}_{type}.bam')
    threads: THREADS_MAXIMUM_CORE
    params: scratch='scratch/STAR-{sample}-{type}'
    run:
        genomedir = os.path.join(TAILSEEKER_DIR, 'refdb', 'level2',
                        CONF['reference_set'][wildcards.sample], 'index.star')
        input = suffix_filter(input)

        if os.path.isdir(params.scratch):
            shutil.rmtree(params.scratch)
        os.makedirs(params.scratch)

        if CONF['performance']['enable_gsnap']:
            unmapped_opts = '--outSAMunmapped None --outReadsUnmapped Fastx '
        else:
            unmapped_opts = '--outSAMunmapped Within KeepPairs --outReadsUnmapped None '

        shell('{STAR_CMD} --runThreadN {threads} --genomeDir {genomedir} \
                --readFilesIn {input[R5.fastq.gz]} {input[R3.fastq.gz]} \
                --readFilesCommand zcat \
                --outFilterType BySJout \
                --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
                --outFilterMismatchNmax 999 \
                --alignIntronMin 15 --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000 --runRNGseed 8809 \
                --outSAMtype BAM Unsorted --sysShell {BASH_CMD} \
                --quantMode TranscriptomeSAM \
                --quantTranscriptomeBan Singleend \
                --outTmpDir {params.scratch}/tmp --outStd BAM_Unsorted \
                --outFileNamePrefix {params.scratch}/ \
                {unmapped_opts} --outMultimapperOrder Random \
                --outSAMmapqUnique 41 > {output.mapped}')
        shell('mv -f "{params.scratch}/Aligned.toTranscriptome.out.bam" \
                     {output.transcriptome}')

        import gzip
        if not CONF['performance']['enable_gsnap']:
            gzip.open(output.unmapped5, 'w')
            gzip.open(output.unmapped3, 'w')
        elif wildcards.type == 'single':
            shell('{BGZIP_CMD} -@ {threads} -c {params.scratch}/Unmapped.out.mate1 \
                    > {output.unmapped5} && rm -f {params.scratch}/Unmapped.out.mate1')
            gzip.open(output.unmapped3, 'w')
        else:
            # STAR sometimes write unmapped FASTQ files in different order for
            # mates. We here reorder them to keep it usable by GSNAP.
            for mateno, outputfile in zip([1, 2], [output.unmapped5, output.unmapped3]):
                shell('''awk '{{ printf("%s\b", $0); n++; if (n%4==0) {{ print; }} }}' \
                            {params.scratch}/Unmapped.out.mate{mateno} | \
                            sort -k1,1 -t'\b' | \
                            awk -F'\b' '{{ \
                                printf("%s\\n%s\\n%s\\n%s\\n", $1, $2, $3, $4); }}' | \
                            {BGZIP_CMD} -@ {threads} -c > {outputfile}''')

                shell('rm -f {params.scratch}/Unmapped.out.mate{mateno}')

        shutil.rmtree(params.scratch)

def inputs_for_GSNAP_alignment(wildcards):
    fastq_filename = 'scratch/unmapped-reads/{sample}-{type}_{{read}}.fastq.gz'.format(
                            sample=wildcards.sample, type=wildcards.type)

    inputs = [fastq_filename.format(read='R5')]
    if wildcards.type == 'paired':
        inputs.append(fastq_filename.format(read='R3'))
    return inputs

rule GSNAP_alignment:
    input: inputs_for_GSNAP_alignment
    output: temp('scratch/alignments/{sample}_GSNAP_{type,[^_.]+}.bam.{part}')
    threads: THREADS_MAXIMUM_CORE
    run:
        genomedir = os.path.join(TAILSEEKER_DIR, 'refdb', 'level2',
                        CONF['reference_set'][wildcards.sample], 'index.gmap')
        partno = '{}/{}'.format(wildcards.part, CONF['performance']['split_gsnap_jobs'])

        shell('{GSNAP_CMD} -D {genomedir} -d genome -A sam -B 4 --gunzip -q {partno} \
                -s {genomedir}/splicesites.iit -m 0.05 -t {threads} \
                {input} | {SAMTOOLS_CMD} view -@ 3 -bS - > {output}')

if CONF['performance']['enable_gsnap']:
    rule merge_alignments:
        input:
            star='scratch/alignments/{sample}_STAR_{type}.bam',
            gsnap=expand('scratch/alignments/{{sample}}_GSNAP_{{type}}.bam.{part}',
                         part=range(CONF['performance']['split_gsnap_jobs'])),
            taginfo=expand('scratch/taginfo/{{sample}}_{tile}.txt.gz', tile=TILES)
        output: temp('scratch/merged-alignments/{sample}_{type,[^_.]+}.bam')
        threads: THREADS_MAXIMUM_CORE
        params: sorttmp='scratch/alignments/{sample}_merge_{type}'
        # samtools 1.3 merge does not respect `-n' option for paired alignments.
        shell: '{SAMTOOLS_CMD} merge -n -u -h {input.star} -@ {threads} - \
                    {input.star} {input.gsnap} | \
                {SAMTOOLS_CMD} sort -n -@ {threads} -T {params.sorttmp} -O sam - | \
                {PYTHON3_CMD} {SCRIPTSDIR}/add-sam-tags-primary.py {input.taginfo} | \
                {SAMTOOLS_CMD} view -@ {threads} -b -o {output} -'
else:
    rule merge_alignments:
        input:
            star='scratch/alignments/{sample}_STAR_{type}.bam',
            taginfo=expand('scratch/taginfo/{{sample}}_{tile}.txt.gz', tile=TILES)
        output: temp('scratch/merged-alignments/{sample}_{type,[^_.]+}.bam')
        threads: THREADS_MAXIMUM_CORE
        params: sorttmp='scratch/alignments/{sample}_merge_{type}'
        shell: '{SAMTOOLS_CMD} sort -n -@ {threads} -T {params.sorttmp} -O sam {input.star} | \
                {PYTHON3_CMD} {SCRIPTSDIR}/add-sam-tags-primary.py {input.taginfo} | \
                {SAMTOOLS_CMD} view -@ {threads} -b -o {output} -'


# ---
# Refinement of non-templated 3'-end additions with reference genome.
# ---

rule tabix_index_taginfo:
    input: 'taginfo/{sample}.txt.gz'
    output: temp('taginfo/{sample}.txt.gz.tbi')
    shell: '{TABIX_CMD} -0 -b 2 -e 2 -s 1 {input}'

TARGETS.extend(expand('refined-taginfo/{sample}.{type}.txt.gz',
                      sample=EXP_SAMPLES, type=['mapped', 'all']))
rule reevaluate_tails:
    input:
        taginfo='taginfo/{sample}.txt.gz',
        taginfo_index='taginfo/{sample}.txt.gz.tbi',
        alignment='scratch/merged-alignments/{sample}_paired.bam'
    output: temp('refined-taginfo/{sample}.txt.pre.gz')
    threads: THREADS_MAXIMUM_CORE
    run:
        genomedir = os.path.join(TAILSEEKER_DIR, 'refdb', 'level2',
                                 CONF['reference_set'][wildcards.sample])
        analytic_options = """\
            --max-fragment-size {modopt[maximum_fragment_size]} \
            --max-polya-reevaluation {modopt[polya_reevaluation_limit]} \
            --terminal-mismatch-recheck {modopt[terminal_mismatch_refine]} \
            --polya-reev-contamination-prob {modopt[polya_reevaluation_contamination_prob]}\
            --contaminated-polya-rescue-threshold \
                {modopt[polya_reevaluation_terminal_a_rescue]}""".format(
                modopt=CONF['modification_refinement'])

        shell('{SCRIPTSDIR}/refine-modifications.py \
                --parallel {threads} \
                --taginfo {input.taginfo} --alignment {input.alignment} \
                --reference-seq {genomedir}/genome.fa {analytic_options} | \
                {BGZIP_CMD} -c > {output}')

rule generate_short_polya_list:
    input: 'refined-taginfo/{sample}.txt.pre.gz'
    output: temp('scratch/short-polya-list/{sample}.txt')
    run:
        import pandas as pd
        from tailseeker import tabledefs

        tbl = pd.read_table(input[0], **tabledefs.refined_taginfo)
        length_lo, length_hi = CONF['modification_refinement']['alignable_polya_range']
        patags = tbl[(tbl['polyA'] >= length_lo) & (tbl['polyA'] <= length_hi)].copy()
        tailids = patags.apply(lambda row: tabledefs.seqid_format.format(r=row), axis=1)
        tailids.to_csv(output[0], header=False, index=False)

rule extract_short_polya_tag_alignments:
    input:
        idlist='scratch/short-polya-list/{sample}.txt',
        alignments='scratch/merged-alignments/{sample}_paired.bam'
    output: temp('scratch/polya-sites/indiv-{sample}.txt')
    threads: 4
    run:
        genomedir = os.path.join(TAILSEEKER_DIR, 'refdb', 'level2',
                                 CONF['reference_set'][wildcards.sample])

        shell('{SAMTOOLS_CMD} view -h -f 128 -F 4 {input.alignments} | \
               {SCRIPTSDIR}/sam-whitelist-filter.py {input.idlist} | \
               {SAMTOOLS_CMD} view -bS - | \
               {BEDTOOLS_CMD} bamtobed -i - | \
               sed -e \'s,-$,x,g\' -e \'s,+$,-,g\' -e \'s,x$,+,g\' | \
               {BEDTOOLS_CMD} flank -s -l 0 -r 1 -i - -g {genomedir}/chrom-sizes | \
               awk -F\'\t\' \'($2 < $3) {{ print $0; }}\' | \
               {BEDTOOLS_CMD} slop -s -l 1 -r -1 -g {genomedir}/chrom-sizes | \
               {BEDTOOLS_CMD} sort -i - > {output}')

def inputs_for_merge_polya_sites_list(wildcards):
    return ['scratch/polya-sites/indiv-{sample}.txt'.format(sample=s)
            for s in EXPERIMENT_GROUPS[wildcards.group]]

rule merge_polya_sites_list:
    input: inputs_for_merge_polya_sites_list
    output: temp('scratch/polya-sites/group-{group}.txt')
    run:
        sample1 = EXPERIMENT_GROUPS[wildcards.group][0]
        genomedir = os.path.join(TAILSEEKER_DIR, 'refdb', 'level2',
                                 CONF['reference_set'][sample1])
        polya_site_window = CONF['modification_refinement']['polya_site_flank']

        shell('cat {input} | {BEDTOOLS_CMD} sort -i - | \
               {BEDTOOLS_CMD} merge -s -i - -c 4,5,6 -o first | \
               {BEDTOOLS_CMD} slop -s -l {polya_site_window} \
                    -r {polya_site_window} -i - -g {genomedir}/chrom-sizes | \
               {BEDTOOLS_CMD} merge -i - -s -c 4,5,6 -o first > {output}')

def inputs_for_apply_short_polya_filter(wildcards):
    return 'scratch/polya-sites/group-{group}.txt'.format(
                group=CONF['experiment_groups'][wildcards.sample])

rule apply_short_polya_filter:
    input:
        alignment='scratch/merged-alignments/{sample}_paired.bam',
        taginfo='refined-taginfo/{sample}.txt.pre.gz',
        pasitelist=inputs_for_apply_short_polya_filter
    output: 'refined-taginfo/{sample}.all.txt.gz'
    params: tmplist=SCRATCHDIR+'/apply_short_polya_filter-{sample}.txt'
    run:
        genomedir = os.path.join(TAILSEEKER_DIR, 'refdb', 'level2',
                                 CONF['reference_set'][wildcards.sample])

        shell("{SAMTOOLS_CMD} view -b -f 128 -F 4 {input.alignment} | \
               {BEDTOOLS_CMD} bamtobed | \
               sed -e 's,-$,x,g' -e 's,+$,-,g' -e 's,x$,+,g' | \
               {BEDTOOLS_CMD} flank -s -l 0 -r 1 -i - -g {genomedir}/chrom-sizes | \
               awk -F '\t' '($2 < $3) {{ print $0; }}' | \
               {BEDTOOLS_CMD} slop -s -l 1 -r -1 -g {genomedir}/chrom-sizes | \
               {BEDTOOLS_CMD} intersect -nonamecheck -u -s -a - -b {input.pasitelist} | \
               cut -f4 | cut -d: -f1-2 | uniq > {params.tmplist}")

        import pandas as pd
        import numpy as np
        import subprocess as sp
        import io
        from tailseeker import tabledefs

        shorttailsid = pd.read_table(params.tmplist, sep=':', names=['tile', 'cluster'],
                                     dtype={'tile': str, 'cluster': np.uint32})
        shorttailsid['additional_flags'] = 1.0

        taginfo = pd.read_table(input.taginfo, **tabledefs.refined_taginfo)
        merged = pd.merge(taginfo, shorttailsid, how='left', left_on=('tile', 'cluster'),
                          right_on=('tile', 'cluster'))
        merged['pflags'] += (tabledefs.PAFLAG_LIKELY_HAVE_INTACT_END *
                             merged['additional_flags'].notnull())

        with sp.Popen([BGZIP_CMD, '-c'], stdin=sp.PIPE, stdout=open(output[0], 'wb')) as proc:
            with io.TextIOWrapper(proc.stdin) as stdin_text_stream:
                merged.iloc[:, :len(taginfo.columns)].to_csv(stdin_text_stream,
                                                        sep='\t', index=False, header=False)

        os.unlink(params.tmplist)


TARGETS.append('stats/polya-length-distributions-L2.csv')
rule generate_polya_length_distribution_stats_level2:
    input: expand('refined-taginfo/{sample}.mapped.txt.gz', sample=sorted(EXP_SAMPLES))
    output: 'stats/polya-length-distributions-L2.csv'
    params:
        samplenames=sorted(EXP_SAMPLES),
        badflagmask=CONF['qcstats']['bad_flags_filter'],
        refined=True, maxpalength=CONF['read_cycles']['R3'][1]
    script: SCRIPTSDIR + '/stats-polya-len-dists.py'


# ---
# Duplicated tag elimination by approximate matching using R5 mapped positions.
# ---

rule sort_alignments:
    input:
        bam='scratch/merged-alignments/{sample}_{type}.bam',
        taginfo='refined-taginfo/{sample}.all.txt.gz'
    output: temp('scratch/sorted-alignments/{sample}_{type,[^_.]+}.bam')
    threads: THREADS_MAXIMUM_CORE
    params: sorttmp='scratch/merged-alignments/{sample}_{type}'
    shell: '{SAMTOOLS_CMD} view -h {input.bam} | \
            {PYTHON3_CMD} {SCRIPTSDIR}/add-sam-tags-refined.py {input.taginfo} | \
            {SAMTOOLS_CMD} sort -@ {threads} -T {params.sorttmp} -o {output} -'

rule index_sorted_alignments:
    input: 'scratch/sorted-alignments/{name}.bam'
    output: temp('scratch/sorted-alignments/{name}.bam.bai')
    shell: '{SAMTOOLS_CMD} index -b {input} {output}'

rule find_approximate_duplicates:
    input:
        bam='scratch/sorted-alignments/{sample}_single.bam',
        bamidx='scratch/sorted-alignments/{sample}_single.bam.bai',
    output: temp('scratch/approx-duplicates/{sample}.txt')
    threads: 6
    run:
        dedupopts = CONF['approximate_duplicate_elimination']
        shell('{BINDIR}/tailseq-dedup-approx {input.bam} \
                {dedupopts[mapped_position_tolerance]} \
                {dedupopts[umi_edit_dist_tolenrance]} {threads} | \
               sort -k3,3 | uniq -f 2 > {output}')

rule filter_approximate_duplicates:
    input:
        bam='scratch/sorted-alignments/{sample}_{type}.bam',
        dupinfo='scratch/approx-duplicates/{sample}.txt'
    output: 'alignments/{sample}_{type,[^_.]+}.bam'
    threads: 3
    shell: '{PYTHON3_CMD} {SCRIPTSDIR}/filter-approximate-duplicates.py \
                --bam {input.bam} --duplicates {input.dupinfo} | \
            {SAMTOOLS_CMD} view -b -@ {threads} -o {output} -'

rule update_refined_taginfo_for_mapped:
    input:
        taginfo='refined-taginfo/{sample}.all.txt.gz',
        dupinfo='scratch/approx-duplicates/{sample}.txt'
    output: 'refined-taginfo/{sample}.mapped.txt.gz'
    run:
        import pandas as pd
        import numpy as np
        import subprocess as sp
        import io
        from tailseeker import tabledefs

        dupinfo = pd.read_table(input.dupinfo, names=['polyA', 'clones', 'readid'],
                                dtype={'polyA': np.int32, 'clones': np.uint32,
                                       'readid': str})
        dupinfo['tile'] = dupinfo['readid'].apply(lambda x: x.split(':')[0])
        dupinfo['cluster'] = dupinfo['readid'].apply(lambda x: int(x.split(':')[1]))

        taginfo = pd.read_table(input.taginfo, **tabledefs.refined_taginfo)
        merged = pd.merge(taginfo, dupinfo, how='right', left_on=('tile', 'cluster'),
                          right_on=('tile', 'cluster'), suffixes=['_orig', ''])

        with sp.Popen([BGZIP_CMD, '-c'], stdin=sp.PIPE, stdout=open(output[0], 'wb')) as proc:
            with io.TextIOWrapper(proc.stdin) as stdin_text_stream:
                merged[taginfo.columns].to_csv(stdin_text_stream, sep='\t',
                                               index=False, header=False)

rule index_alignments:
    input: 'alignments/{name}.bam'
    output: 'alignments/{name}.bam.bai'
    shell: '{SAMTOOLS_CMD} index -b {input} {output}'

rule make_read5_position_dists:
    input:
        bam='scratch/tr-alignments/{sample}_single.bam',
        transcript_sizes=lambda wc: os.path.join(TAILSEEKER_DIR, 'refdb', 'level2',
                                        CONF['reference_set'][wc.sample], 'transcript-sizes')
    output: temp('scratch/read5-pos-dists/{sample}.csv')
    script: SCRIPTSDIR + '/make-read5-positions-dist.py'


# ---
# Quality check and stats
# ---

TARGETS.append('stats/read5-position-dists.csv')

rule merge_read5_position_dists:
    input: expand('scratch/read5-pos-dists/{sample}.csv', sample=EXP_SAMPLES)
    output: 'stats/read5-position-dists.csv'
    run:
        import pandas, numpy
        (pandas.DataFrame({
            name: pandas.read_csv(infile, names=['pos', 'count'], index_col=0)['count']
            for name, infile in zip(EXP_SAMPLES, input)})
         .fillna(0).astype(numpy.int64).to_csv(output[0]))


# ---
# Associations to genes and gene-level statistics
# ---

TARGETS.extend(expand('associations/{sample}.txt.gz', sample=EXP_SAMPLES))

rule associate_tags_to_genes:
    input: 'alignments/{sample}_single.bam'
    output: temp('scratch/associations/{sample}.txt')
    threads: 4
    run:
        genomedir = os.path.join(TAILSEEKER_DIR, 'refdb', 'level2',
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
                "sed -e \'s,^\\([^:]*\\):0*\\([0-9][0-9]*\\)\\(:[^\t]*\\)*\t,\\1\t\\2\t,g\'" | \
            {BGZIP_CMD} -@ {threads} -c > {output}'


TARGETS.extend(expand('tagcounts/{sample}-{ambigtype}-{modtype}.msgpack.xz',
                      sample=EXP_SAMPLES, ambigtype=['single', 'multi'],
                      modtype=['U', 'C', 'G']))

rule make_gene_level_counts:
    input:
        taginfo='refined-taginfo/{sample}.mapped.txt.gz',
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


TARGETS.extend(expand('stats/genelevelstats-{sample}-{ambigtype}.csv',
                      sample=EXP_SAMPLES, ambigtype=['single', 'multi']))

rule make_gene_level_statistics:
    input:
        U='tagcounts/{sample}-{ambigtype}-U.msgpack.xz',
        G='tagcounts/{sample}-{ambigtype}-G.msgpack.xz',
        C='tagcounts/{sample}-{ambigtype}-C.msgpack.xz'
    output: limit('stats/genelevelstats-{sample}-{{ambigtype}}.csv', sample=EXP_SAMPLES)
    params: confidence_interval_span=CONF['gene_level_stats']['polyA_len_confidence_interval']
    script: SCRIPTSDIR + '/stats-gene-level-tailing.py'

# ex: syntax=snakemake
