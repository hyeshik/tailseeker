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

ALL_GENOMES_USED = sorted(set(CONF['reference_set'].values()))
SAMPLES_BY_GENOME = {
    genome: [smp for smp in EXP_SAMPLES if CONF['reference_set'][smp] == genome]
    for genome in ALL_GENOMES_USED
}

TARGETS.extend(expand('stats/genelevelstats-{genome}-{ambigtype}.{ext}',
                      genome=ALL_GENOMES_USED, ambigtype=['single', 'multi'],
                      ext=['csv', 'feather', 'xlsx']))

rule merge_gene_level_stats:
    input:
        lambda wc: expand('stats/genelevelstats-{sample}-{ambigtype}.csv',
                          sample=SAMPLES_BY_GENOME[wc.genome], ambigtype=[wc.ambigtype])
    output:
        csv=limit('stats/genelevelstats-{genome}-{{ambigtype}}.csv', genome=ALL_GENOMES_USED),
        feather='stats/genelevelstats-{genome}-{ambigtype}.feather',
        excel='stats/genelevelstats-{genome}-{ambigtype}.xlsx'
    params:
        genomedir=REFDBDIR + '/level3/{genome}',
        samples_by_genome=SAMPLES_BY_GENOME
    script: SCRIPTSDIR + '/merge-gene-level-stats.py'

# ex: syntax=snakemake
