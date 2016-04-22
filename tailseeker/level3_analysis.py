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
    output: 'associations/{sample}.txt.gz'
    threads: 6
    run:
        genomedir = os.path.join(TAILSEEKER_DIR, 'refdb', 'level3',
                                 CONF['reference_set'][wildcards.sample])
        required_mapq = CONF['gene_level_stats']['required_mapping_quality']

        shell('{SAMTOOLS_CMD} view -b -F 260 -q {required_mapq} {input} | \
               {BEDTOOLS_CMD} intersect -nonamecheck -split -wa -wb -bed \
                -a - -b {genomedir}/exons.gtf.gz -s | \
               cut -f4,21 | \
               {PARALLEL_CMD} -j {threads} --pipe --block 10M \
                    "sed -e \'s,gene_id \\"\\([^\\"]*\\)\\".*$,\\1,g\' \
                        -e \'s,^\\([^:]*\\):0*\\([0-9][0-9]*\\):[^\t]*\t,\\1\t\\2\t,g\'" | \
               sort -k1,1 -k2,2n -k3,3 | uniq | \
               {BGZIP_CMD} -c > {output}')

# ex: syntax=snakemake
