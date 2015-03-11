#!/usr/bin/env python3
#
# Copyright (c) 2013-2015 Institute for Basic Science
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

# TODO: move these to global configuration
THREADS_MAXIMUM_CORE = 32
AYB_BINDIR = '/atp/hyeshik/gh/AYB2/bin'
HTSLIB_BINDIR = '/usr/local/bin'
#=======================

TARGETS = []

include: 'tailor/workflows/snakesupport.py'

from workflows import sequencers

# Variable settings
TILES = sequencers.get_tiles(CONF)

EXP_SAMPLES = sorted(CONF['experimental_samples'].keys())
SPIKEIN_SAMPLES = sorted(CONF['spikein_samples'].keys())
ALL_SAMPLES = sorted(EXP_SAMPLES + SPIKEIN_SAMPLES)

INSERT_READS = sorted(readname for readname in CONF['read_cycles'] if readname[1] != 'i')
INDEX_READS = sorted(readname for readname in CONF['read_cycles'] if readname[1] == 'i')
ALL_READS = INSERT_READS + INDEX_READS

FIRST_CYCLE = min(f for f, l in CONF['read_cycles'].values())
LAST_CYCLE = max(l for f, l in CONF['read_cycles'].values())
NUM_CYCLES = LAST_CYCLE - FIRST_CYCLE + 1

PHIX_ID_REF = ['R5', 6, 40] # identify PhiX reads using 40 bases from the 6th cycle.

# Variable validations
if len(INDEX_READS) != 1:
    raise ValueError("Multi-indexing is not supported yet.")

if FIRST_CYCLE != 1:
    raise ValueError("The pipeline assumes that one of the reads starts from the first cycle.")


localrules: all

rule all:
    input: lambda wc: TARGETS

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
        import shutil

        tileinfo = TILES[wildcards.tile]
        readname = wildcards.read
        first_cycle, last_cycle = CONF['read_cycles'][readname]
        read_length = last_cycle - first_cycle + 1

        tempdir = make_scratch_dir('aybcalls/{}_{}'.format(readname, tileinfo['id']))
        reads_format = (('' if first_cycle == 1 else 'I{}'.format(first_cycle - 1)) +
                        'R{}'.format(read_length))

        shell('{AYB_BINDIR}/AYB -p {threads} -o {tempdir} -i {tileinfo[datadir]} \
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

rule demultiplex_signals:
    """
    Extract signals and basecalls from .CIF, .BCL, and FASTQ files. The sequence and signals
    are sorted into files in the in-house .SQI format after demultiplexing to separate
    reads by their indices (barcodes).
    """
    input: determine_inputs_demultiplex_signals
    output:
        regular = map(temp, expand('scratch/demux-sqi/{sample}_{{tile}}.sqi.gz',
                                   sample=ALL_SAMPLES)),
        unknown = temp('scratch/demux-sqi/Unknown_{tile}.sqi.gz'),
        phixcontrol = temp('scratch/demux-sqi/PhiX_{tile}.sqi.gz')
    threads: min(len(EXP_SAMPLES) + 2, 8)
    run:
        tileinfo = TILES[wildcards.tile]
        index_read_start, index_read_end = CONF['read_cycles'][INDEX_READS[0]]
        index_read_length = index_read_end - index_read_start + 1

        output_filename = 'scratch/demux-sqi/XX_{tile}.sqi.gz'.format(tile=wildcards.tile)

        phix_cycle_base = CONF['read_cycles'][PHIX_ID_REF[0]][0]
        phix_cycle_start = phix_cycle_base + PHIX_ID_REF[1] - 1
        phix_cycle_length = PHIX_ID_REF[2]

        options = [
            '--data-dir', os.path.join(tileinfo['datadir'], 'Data', 'Intensities'),
            '--run-id', tileinfo['laneid'], '--lane', tileinfo['lane'],
            '--tile', tileinfo['tile'], '--ncycles', NUM_CYCLES,
            '--signal-scale', sequencers.get_signalscale(tileinfo['type']),
            '--barcode-start', index_read_start, '--barcode-length', index_read_length,
            '--writer-command', '{HTSLIB_BINDIR}/bgzip -c > ' + output_filename,
            '--filter-control', 'PhiX,{},{}'.format(phix_cycle_start, phix_cycle_length),
        ]

        altcalls = determine_inputs_demultiplex_signals(wildcards, as_dict=True)
        for readid, inputfile in altcalls.items():
            options.extend([
                '--alternative-call',
                '{},{}'.format(inputfile, CONF['read_cycles'][readid][0])])

        for sample in ALL_SAMPLES:
            options.extend([
                '--sample',
                '{name},{index},{mismatches},{delimiter},{delimpos}'.format(
                    name=sample, index=CONF.get_sample_index(sample),
                    mismatches=CONF['maximum_index_mismatches'][sample],
                    delimiter=CONF['delimiter'][sample][1],
                    delimpos=CONF['delimiter'][sample][0])])

        shell('{BINDIR}/tailseq-retrieve-signals ' +
              ' '.join('"{}"'.format(opt) for opt in options))


rule merge_sqi:
    input: expand('scratch/demux-sqi/{{sample}}_{tile}.sqi.gz', tile=TILES)
    output: 'scratch/merged-sqi/{sample}.sqi.gz'
    run:
        input = sorted(input) # to make tabix happy.
        shell('{SCRIPTSDIR}/bgzf-merge.py --output {output} {input}')

TARGETS.extend(expand('scratch/merged-sqi/{sample}.sqi.gz', sample=ALL_SAMPLES))
TARGETS.extend(expand('scratch/demux-sqi/{sample}_{tile}.sqi.gz', sample=['Unknown', 'PhiX'],
                      tile=TILES))

# ex: syntax=snakemake
