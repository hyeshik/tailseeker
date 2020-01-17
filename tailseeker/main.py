#!/usr/bin/env python3
#
# Copyright (c) 2013-2016 Hyeshik Chang
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

import os
import sys
import shutil

TAILSEEKER_DIR = os.path.abspath(os.path.join(os.path.dirname(workflow.snakefile), '..'))
if TAILSEEKER_DIR not in sys.path:
    sys.path.append(TAILSEEKER_DIR)
TARGETS = []

include: os.path.join(TAILSEEKER_DIR, 'tailseeker', 'snakesupport.py')

from tailseeker import sequencers

# Variable settings
TILES = sequencers.get_tiles(CONF)
for badtile in CONF['input_filtering']['blacklisted_tiles']:
    if badtile in TILES:
        del TILES[badtile]

EXP_SAMPLES = CONF.exp_samples
SPIKEIN_SAMPLES = CONF.spikein_samples
ALL_SAMPLES = sorted(EXP_SAMPLES + SPIKEIN_SAMPLES)

INSERT_READS = sorted(readname for readname in CONF['read_cycles'] if readname[1] != 'i')
INDEX_READS = sorted(readname for readname in CONF['read_cycles'] if readname[1] == 'i')
ALL_READS = INSERT_READS + INDEX_READS

FIRST_CYCLE = min(f for f, l, _ in CONF['read_cycles'].values())
LAST_CYCLE = max(l for f, l, _ in CONF['read_cycles'].values())
NUM_CYCLES = LAST_CYCLE - FIRST_CYCLE + 1

THREADS_MAXIMUM_CORE = CONF['maximum_threads']

INTERMEDIATE_DIRS = [
    'dupfilter', 'polya', 'scores', 'scratch',
    'sequences', 'signalproc', 'learning',
]


# Variable validations
if len(INDEX_READS) != 1:
    raise ValueError("Multi-indexing is not supported yet.")

if FIRST_CYCLE != 1:
    raise ValueError("The pipeline assumes that one of the reads starts from the first cycle.")


def sorted_spikein_first(samples):
    spikein_lengths = CONF['spikein_lengths'] if 'spikein_lengths' in CONF else {}

    names_with_design_length = [
        (spikein_lengths.get(sample, inf), sample)
        for sample in samples]

    return [name for _, name in sorted(names_with_design_length)]


localrules: all

rule all:
    input: lambda wc: TARGETS


def clear_generated_files(include_targets=False):
    cleared_files = 0
    for pattern in all_intermediate_files:
        from snakemake.utils import listfiles

        for filename, wildcards in listfiles(pattern):
            if os.path.isfile(filename):
                os.unlink(filename)
                cleared_files += 1
    else:
        if cleared_files > 0:
            print('{} intermediate file(s) were removed.'.format(cleared_files))

    for dir in INTERMEDIATE_DIRS:
        if os.path.islink(dir):
            realpath = os.readlink(dir)
            print('Cleaning a symbolic link:', dir)
            shutil.rmtree(realpath, ignore_errors=True)
            os.unlink(dir)
        elif os.path.isdir(dir):
            print('Cleaning an intermediate directory:', dir)
            shutil.rmtree(dir, ignore_errors=True)

    if not include_targets:
        return

    from itertools import groupby

    files_safe_to_delete = {'.AppleDouble', '.DS_Store'}

    for dirname, files in groupby(sorted(TARGETS), lambda x: x.split(os.sep)[0]):
        for f in files:
            if os.path.exists(f):
                os.unlink(f)

        if os.path.isdir(dirname) and not (
                set(os.listdir(dirname)) - files_safe_to_delete):
            shutil.rmtree(dirname)

rule clean:
    run:
        clear_generated_files(include_targets=False)

rule clear:
    run:
        clear_generated_files(include_targets=True)


rule basecall_ayb:
    """
    Runs a third-party basecaller, AYB, to call bases. In most cases, the sequence reads
    determined by AYB for 5' end of inserts are more readily aligned to the genome than
    performing that using Illumina RTA.

    Download my modified tree from https://github.com/hyeshik/AYB2. There is few bug fixes
    that were not incorporated into the original version.
    """
    output: temp('scratch/aybcalls/{read}_{tile}.fastq.gz')
    threads: THREADS_MAXIMUM_CORE
    run:
        tileinfo = TILES[wildcards.tile]
        readname = wildcards.read
        first_cycle, last_cycle, read_no = CONF['read_cycles'][readname]
        read_length = last_cycle - first_cycle + 1

        tempdir = make_scratch_dir('aybcalls/{}_{}'.format(readname, tileinfo['id']))
        reads_format = (('' if first_cycle == 1 else 'I{}'.format(first_cycle - 1)) +
                        'R{}'.format(read_length))

        shell('{AYB_CMD} -p {threads} -o {tempdir} -i {tileinfo[datadir]} \
                -b {reads_format} -f fastq.gz -r L{tileinfo[lane]}T{tileinfo[tile]}')
        shell('mv {tempdir}/s_{tileinfo[lane]}_{tileinfo[tile]}.fastq.gz {output}')
        shutil.rmtree(tempdir)


def determine_inputs_process_signals(wildcards):
    inputs = []

    for readid, program in CONF.get('third_party_basecaller', {}).items():
        if program:
            fastqpath = 'scratch/{program}calls/{readid}_{tile}.fastq.gz'.format(
                            program=program.lower(), readid=readid, tile=wildcards.tile)
            inputs.append(fastqpath)

    return inputs

for readid, program in CONF.get('third_party_basecaller', {}).items():
    if program and program not in CONF.paths:
        logger.error('{} is not configured in this installation. '
                     'Please run setup.sh again.'.format(program))
        raise ValueError('Configuration error: {} not in paths.conf'.format(program))


rule process_signals:
    input: determine_inputs_process_signals
    output:
        sigproc_conf = temp('scratch/sigproc-conf/{tile}.ini'),
        seqqual = map(temp, expand('scratch/seqqual/{sample}_{{tile}}.txt.gz',
                                   sample=ALL_SAMPLES)),
        taginfo = map(temp, expand('scratch/taginfo/{sample}_{{tile}}.txt.gz',
                                   sample=ALL_SAMPLES)),
        signals = map(temp, expand('scratch/signals/{sample}_{{tile}}.sigpack',
                                   sample=ALL_SAMPLES)),
        sigdists = map(temp, expand('scratch/sigdists-r00/{posneg}_{{tile}}.sigdists',
                                    posneg=['pos', 'neg'])),
        demuxstats = temp('scratch/stats/signal-proc-{tile}.csv')
    threads: THREADS_MAXIMUM_CORE
    params:
        tileinfo=TILES, conf=CONF.confdata,
        exp_samples=EXP_SAMPLES, spikein_samples=SPIKEIN_SAMPLES
    run:
        external_script('{PYTHON3_CMD} {SCRIPTSDIR}/generate-signalproc-conf.py')
        shell('{BINDIR}/tailseq-import {output.sigproc_conf}')


TARGETS.append('stats/signal-processing.csv')
rule merge_signal_processing_stats:
    input: expand('scratch/stats/signal-proc-{tile}.csv', tile=sorted(TILES))
    output: 'stats/signal-processing.csv'
    params: tiles=sorted(TILES)
    script: SCRIPTSDIR + '/stats-merge-demultiplexing-counts.py'


rule generate_polyA_finder_init_parameters:
    input:
        sigdist_pos = 'scratch/sigdists-r00/pos_{tile}.sigdists',
        sigdist_neg = 'scratch/sigdists-r00/neg_{tile}.sigdists'
    output: temp('scratch/avgs_r00/{tile}.avg')
    params: CONF=CONF.confdata
    run: external_script('{PYTHON3_CMD} {SCRIPTSDIR}/generate-init-average.py')
        
rule polyA_find_iterator:
    input:
        avg = lambda wc: ('scratch/avgs_r{rnd:02d}/{tile}.avg'
                          .format(rnd=int(wc.round)-1, tile="{tile}")),
        sigpacks = expand('scratch/signals/{sample}_{{tile}}.sigpack', sample=EXP_SAMPLES)
    output:
        avg = 'scratch/avgs_r{round,[^0].|.[^0]}/{tile}.avg'
    threads: THREADS_MAXIMUM_CORE
    params: CONF=CONF.confdata, isWrite=False
    run: external_script('{PYTHON3_CMD} {SCRIPTSDIR}/polya-find-iterator.py')
    
        
rule polyA_taginfo_writer:
    input:
        avg = 'scratch/avgs_r{rnd:02d}/{tile}.avg'
                     .format(rnd=CONF['polyA_ruler']['signal_resampling_rounds'], tile="{tile}"),
        sigpacks = expand('scratch/signals/{sample}_{{tile}}.sigpack', sample=ALL_SAMPLES),
        taginfo = expand('scratch/taginfo/{sample}_{{tile}}.txt.gz', sample=ALL_SAMPLES)
    output:
        taginfo = map(temp, expand('scratch/taginfo_final/{sample}_{{tile}}.txt.gz', sample=ALL_SAMPLES))
    threads: THREADS_MAXIMUM_CORE
    params: CONF=CONF.confdata, isWrite=True, tileinfo=lambda wc: wc.tile
    run: external_script('{PYTHON3_CMD} {SCRIPTSDIR}/polya-find-iterator.py')
        

TARGETS.extend(expand('taginfo/{sample}.txt.gz', sample=ALL_SAMPLES))
rule merge_and_deduplicate_taginfo:
    # Round number for taginfo files are actually based on parameters
    # of one round behind. Thus, use the file for one step later.
    input:
        expand('scratch/taginfo_final/{{sample}}_{tile}.txt.gz', tile=TILES)
    output:
        taginfo='taginfo/{sample}.txt.gz',
        duptrace=temp('scratch/stats/perfdup-traces-{sample}.txt.gz')
    threads: 12
    run:
        sorted_input = sorted(input)
        if wildcards.sample in EXP_SAMPLES:
            shell('zcat {sorted_input} | \
                env BGZIP_OPT="-@ {threads}" sort -t "\t" -k6,6 -k1,1 -k2,2n \
                    --compress-program={BINDIR}/bgzip-wrap --parallel={threads} | \
                {BINDIR}/tailseq-dedup-perfect {output.duptrace} | \
                env BGZIP_OPT="-@ {threads}" sort -t "\t" -k1,1 -k2,2n \
                    --compress-program={BINDIR}/bgzip-wrap --parallel={threads} | \
                {BGZIP_CMD} -@ {threads} -c > {output.taginfo}')
        elif wildcards.sample in SPIKEIN_SAMPLES:
            external_script('{PYTHON3_CMD} {SCRIPTSDIR}/bgzf-merge.py')
            shell('echo -n "" | gzip -c - > {output.duptrace}')


if CONF['analysis_level'] >= 2 and CONF['read_filtering']['contaminant_filtering']:
    temp_primary_fastq = temp
    INTERMEDIATE_DIRS.append('fastq')
else:
    temp_primary_fastq = lazy_clearing
    TARGETS.extend(expand('fastq/{sample}_{read}.fastq.gz', sample=EXP_SAMPLES,
                          read=INSERT_READS))

rule produce_fastq_outputs:
    input:
        taginfo='taginfo/{sample}.txt.gz',
        seqquals=expand('scratch/seqqual/{{sample}}_{tile}.txt.gz', tile=TILES)
    output:
        R5=temp_primary_fastq('fastq/{sample}_R5.fastq.gz'),
        R3=temp_primary_fastq('fastq/{sample}_R3.fastq.gz')
    threads: THREADS_MAXIMUM_CORE
    params: seqqual_filename='scratch/seqqual/{sample}_@tile@.txt.gz'
    run:
        verbosity_opt = '--fastq-id-verbose ' if CONF['analysis_level'] <= 1 else ''
        shell('{BINDIR}/tailseq-writefastq \
                --taginfo {input.taginfo} --seqqual \'{params.seqqual_filename}\' \
                --fastq5 {output.R5} --fastq3 {output.R3} \
                {verbosity_opt} --threads {threads}')


TARGETS.append('stats/polya-length-distributions-L1.csv')
rule generate_polya_length_distribution_stats_level1:
    input: expand('taginfo/{sample}.txt.gz', sample=sorted_spikein_first(ALL_SAMPLES))
    output: 'stats/polya-length-distributions-L1.csv'
    params:
        samplenames=sorted_spikein_first(ALL_SAMPLES),
        badflagmask=CONF['qcstats']['bad_flags_filter'],
        refined=False, maxpalength=CONF['read_cycles']['R3'][1]
    script: SCRIPTSDIR + '/stats-polya-len-dists.py'


TARGETS.append('qcplots/global-polya-length.pdf')
rule plot_global_polya_length_virtual_gel:
    input: 'stats/polya-length-distributions-L1.csv'
    output: 'qcplots/global-polya-length.pdf'
    params:
        controls=','.join(name for length, name
                          in sorted((CONF['spikein_lengths'][s], s) for s in SPIKEIN_SAMPLES)
                          if length > 0),
        kde_bandwidth=CONF['qcstats']['virtual_gel_kde_bandwidth'],
        minimum_polya_length=CONF['qcstats']['virtual_gel_minimum_polya']
    run:
        samples = '--samples ' + ','.join(EXP_SAMPLES) if EXP_SAMPLES else ''
        shell('{PYTHON3_CMD} {SCRIPTSDIR}/plot-virtual-gel.py \
                    --tagcounts {input} --controls {params.controls} \
                    --minimum-polya-length {params.minimum_polya_length} \
                    --kde-bandwidth {params.kde_bandwidth} \
                    {samples} --output-plot {output}')


TARGETS.append('qcplots/global-polya-length-histogram-L1.pdf')
rule plot_global_polya_length_histogram:
    input: 'stats/polya-length-distributions-L1.csv'
    output: 'qcplots/global-polya-length-histogram-L1.pdf'
    params:
        exclude=SPIKEIN_SAMPLES,
        kde_bandwidth=CONF['qcstats']['histogram_kde_bandwidth'],
        minimum_polya_length=CONF['qcstats']['histogram_minimum_polya'],
        x_transform_factor=CONF['qcstats']['histogram_xscale_factor']
    script: SCRIPTSDIR + '/plot-polya-len-histogram.py'


if 'spikein_lengths' in CONF and len(CONF['spikein_lengths']) > 0:
    TARGETS.extend(['stats/control-length-accuracy.csv',
                    'qcplots/control-length-accuracy.pdf'])

    rule generate_accuracy_stats:
        input: 'stats/polya-length-distributions-L1.csv'
        output:
            statsout='stats/control-length-accuracy.csv',
            plotout='qcplots/control-length-accuracy.pdf'
        run:
            controlsamples = ' '.join('{}:{}'.format(CONF['spikein_lengths'][s], s)
                                      for s in SPIKEIN_SAMPLES)

            shell('{PYTHON3_CMD} {SCRIPTSDIR}/plot-polya-calls-accuracy.py \
                        --control {controlsamples} --output-plot {output.plotout} \
                        --output-stats {output.statsout} \
                        --input {input}')


if CONF['analysis_level'] >= 2:
    include: os.path.join(TAILSEEKER_DIR, 'tailseeker', 'level2_analysis.py')

if CONF['analysis_level'] >= 3:
    include: os.path.join(TAILSEEKER_DIR, 'tailseeker', 'level3_analysis.py')

# ex: syntax=snakemake
