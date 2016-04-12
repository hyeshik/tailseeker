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
import shutil

def get_topdir():
    if os.path.islink('Snakefile'):
        tailseekerpkgdir = os.path.dirname(os.readlink('Snakefile'))
        return os.path.abspath(os.path.dirname(tailseekerpkgdir))
    elif 'TAILSEEKER_DIR' in os.environ:
        return os.environ['TAILSEEKER_DIR']
    else:
        raise ValueError("You need to set an environment variable, TAILSEEKER_DIR.")

TAILSEEKER_DIR = get_topdir()
TARGETS = []

include: os.path.join(TAILSEEKER_DIR, 'tailseeker', 'snakesupport.py')

from tailseeker import sequencers

# Variable settings
TILES = sequencers.get_tiles(CONF)

EXP_SAMPLES = sorted(CONF['experimental_samples'].keys())
SPIKEIN_SAMPLES = sorted(CONF['spikein_samples'].keys())
ALL_SAMPLES = sorted(EXP_SAMPLES + SPIKEIN_SAMPLES)

INSERT_READS = sorted(readname for readname in CONF['read_cycles'] if readname[1] != 'i')
INDEX_READS = sorted(readname for readname in CONF['read_cycles'] if readname[1] == 'i')
ALL_READS = INSERT_READS + INDEX_READS

FIRST_CYCLE = min(f for f, l, _ in CONF['read_cycles'].values())
LAST_CYCLE = max(l for f, l, _ in CONF['read_cycles'].values())
NUM_CYCLES = LAST_CYCLE - FIRST_CYCLE + 1

PHIX_ID_REF = ['R5', 6, 40] # identify PhiX reads using 40 bases from the 6th cycle.

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
    spikein_lengths = CONF['spikein_lengths']

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


def determine_inputs_demultiplex_signals(wildcards, as_dict=False):
    inputs = {}

    for readid, program in CONF.get('third_party_basecaller', {}).items():
        fastqpath = 'scratch/{program}calls/{readid}_{tile}.fastq.gz'.format(
                        program=program.lower(), readid=readid, tile=wildcards.tile)
        inputs[readid] = fastqpath

    return inputs if as_dict else list(inputs.values())

rule process_signals:
    input: determine_inputs_demultiplex_signals
    output:
        sigproc_conf = 'scratch/sigproc-conf/{tile}.ini',
        fastq = map(temp, expand('scratch/fastq/{sample}_{{tile}}_{read}.fastq.gz',
                                 sample=ALL_SAMPLES + ['Unknown', 'PhiX'],
                                 read=INSERT_READS)),
        taginfo = map(temp, expand('scratch/taginfo/{sample}_{{tile}}.txt.gz',
                                   sample=ALL_SAMPLES + ['Unknown', 'PhiX'])),
        demuxstats = temp('scratch/stats/signal-proc-{tile}.csv')
    threads: THREADS_MAXIMUM_CORE
    params:
        tileinfo=TILES, conf=CONF.confdata,
        exp_samples=EXP_SAMPLES, spikein_samples=SPIKEIN_SAMPLES
    run:
        external_script('{PYTHON3_CMD} {SCRIPTSDIR}/generate-signalproc-conf.py')
        shell('{BINDIR}/tailseq-sigproc {output.sigproc_conf}')


TARGETS.extend(expand('scratch/taginfo-sorted/{sample}.txt.gz', sample=EXP_SAMPLES))
rule merge_and_sort_taginfo:
    input: expand('scratch/taginfo/{{sample}}_{tile}.txt.gz', tile=TILES)
    output: 'scratch/taginfo-sorted/{sample}.txt.gz'
    threads: 12
    shell: 'cat {input} | {BGZIP_CMD} -cd -@ {threads} | \
            env BGZIP_OPT="-@ {threads}" sort -t "\t" -k6,6 -k1,1 -k2,2n \
                --compress-program={BINDIR}/bgzip-wrap --parallel={threads} | \
            {BGZIP_CMD} -@ {threads} -c > {output}'

# ex: syntax=snakemake
