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
# Sequence alignment to the genome.
# ---

TARGETS.extend(expand('alignments/{sample}_{type}.bam',
                      sample=EXP_SAMPLES, type=['single', 'paired']))

def inputs_for_STAR_alignment(wildcards):
    fastq_filename = 'fastq/{sample}_{{read}}.fastq.gz'.format(sample=wildcards.sample)

    inputs = [fastq_filename.format(read='R5')]
    if wildcards.type == 'paired':
        inputs.append(fastq_filename.format(read='R3'))
    return inputs

rule STAR_alignment:
    input: inputs_for_STAR_alignment
    output:
        mapped='scratch/alignments/{sample}_STAR_{type,[^_.]+}.bam',
        unmapped='scratch/STAR-{sample}-{type}/Unmapped.out.mate1'
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
                --outTmpDir {params.scratch}/tmp --outStd BAM_Unsorted \
                --outFileNamePrefix {params.scratch}/ \
                {unmapped_opts} --outMultimapperOrder Random \
                --outSAMmapqUnique 41 > {output.mapped}')

        if not CONF['performance']['enable_gsnap']:
            open(output.unmapped, 'w')


rule GSNAP_alignment:
    input: 'scratch/STAR-{sample}-{type}/Unmapped.out.mate1'
    output: 'scratch/alignments/{sample}_GSNAP_{type,[^_.]+}.bam.{part}'
    threads: THREADS_MAXIMUM_CORE
    run:
        genomedir = os.path.join(TAILSEEKER_DIR, 'refdb', 'level2',
                        CONF['reference_set'][wildcards.sample], 'index.gmap')
        if wildcards.type == 'paired':
            input.append(re.sub('mate1$', 'mate2', str(input)))
        partno = '{}/{}'.format(wildcards.part, CONF['performance']['split_gsnap_jobs'])

        shell('{GSNAP_CMD} -D {genomedir} -d genome -A sam -B 4 -q {partno} \
                -s {genomedir}/splicesites.iit -m 0.05 -t {threads} \
                {input} | {SAMTOOLS_CMD} view -@ 3 -bS - > {output}')

if CONF['performance']['enable_gsnap']:
    rule merge_alignments:
        input:
            star='scratch/alignments/{sample}_STAR_{type}.bam',
            gsnap=expand('scratch/alignments/{{sample}}_GSNAP_{{type}}.bam.{part}',
                         part=range(CONF['performance']['split_gsnap_jobs']))
        output: 'alignments/{sample}_{type,[^_.]+}.bam'
        threads: THREADS_MAXIMUM_CORE
        # samtools 1.3 merge does not respect `-n' option for paired alignments.
        shell: '{SAMTOOLS_CMD} merge -n -u -h {input.star} -@ {threads} - {input} | \
                {SAMTOOLS_CMD} sort -n -@ {threads} -O sam - | \
                {PARALLEL_CMD} -j {threads} --keep-order \
                    --pipe "sed -e \'s,^\\([^@][^:]*\\)\\(:.*\\),\\1\\2\tRG:Z:\\1,g\'" | \
                {SAMTOOLS_CMD} view -@ {threads} -b -o {output} -'
else:
    rule merge_alignments:
        input: 'scratch/alignments/{sample}_STAR_{type}.bam'
        output: 'alignments/{sample}_{type,[^_.]+}.bam'
        threads: THREADS_MAXIMUM_CORE
        shell: '{SAMTOOLS_CMD} sort -n -@ {threads} -O sam {input} | \
                {PARALLEL_CMD} -j {threads} --keep-order \
                    --pipe "sed -e \'s,^\\([^@][^:]*\\)\\(:.*\\),\\1\\2\tRG:Z:\\1,g\'" | \
                {SAMTOOLS_CMD} view -@ {threads} -b -o {output} -'


# ---
# Refinement of non-templated 3'-end additions with reference genome.
# ---

rule tabix_index_taginfo:
    input: 'taginfo/{sample}.txt.gz'
    output: temp('taginfo/{sample}.txt.gz.tbi')
    shell: '{TABIX_CMD} -0 -b 2 -e 2 -s 1 {input}'

TARGETS.extend(expand('refined-taginfo/{sample}.txt.gz', sample=EXP_SAMPLES))
rule reevaluate_tails:
    input:
        taginfo='taginfo/{sample}.txt.gz',
        taginfo_index='taginfo/{sample}.txt.gz.tbi',
        alignment='alignments/{sample}_paired.bam'
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
                bgzip -c > {output}')

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
        alignments='alignments/{sample}_paired.bam'
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
        alignment='alignments/{sample}_paired.bam',
        taginfo='refined-taginfo/{sample}.txt.pre.gz',
        pasitelist=inputs_for_apply_short_polya_filter
    output: 'refined-taginfo/{sample}.txt.gz'
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

# ex: syntax=snakemake
