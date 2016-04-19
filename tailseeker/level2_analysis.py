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

        shell('rm -rf {params.scratch}; mkdir -p {params.scratch} && \
               {STAR_CMD} --runThreadN {threads} --genomeDir {genomedir} \
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
                --outReadsUnmapped Fastx --outMultimapperOrder Random \
                --outSAMunmapped None --outSAMmapqUnique 41 > {output}')

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
        shell: '{SAMTOOLS_CMD} merge -u -h {input.star} -@ {threads} - {input} | \
                {SAMTOOLS_CMD} sort -n -@ {threads} -o {output}'
else:
    rule merge_alignments:
        input: 'scratch/alignments/{sample}_STAR_{type}.bam'
        output: 'alignments/{sample}_{type,[^_.]+}.bam'
        threads: THREADS_MAXIMUM_CORE
        shell: '{SAMTOOLS_CMD} sort -n -@ {threads} -o {output} {input}'

# ex: syntax=snakemake
