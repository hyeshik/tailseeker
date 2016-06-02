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

from tailseeker.powersnake import *
import os
import os.path


def generate_source_section(outf):
    tileinfo = params.tileinfo[wildcards.tile]
    read3_rn = params.conf['read_cycles']['R3'][2]
    read3_colormatrix = (
        '{ti[intensitiesdir]}/BaseCalls/Matrix/s_{ti[lane]}_{read3_rn}_'
        '{ti[tile]}_matrix.txt'.format(ti=tileinfo, read3_rn=read3_rn))

    print("""\
[source]
data-dir = {tileinfo[intensitiesdir]}
laneid = {tileinfo[laneid]}
lane = {tileinfo[lane]}
tile = {tileinfo[tile]}
threep-colormatrix = {R3_colormatrix}
""".format(
        tileinfo=tileinfo, R3_colormatrix=read3_colormatrix), file=outf)


def generate_read_format_section(outf):
    read_cycles = params.conf['read_cycles']
    read_length = lambda start, end, rn: end - start + 1
    computed_values = {
        'total_cycles': max(end for _, end, _ in read_cycles.values()),
        'fivep_length': read_length(*read_cycles['R5']),
        'threep_length': read_length(*read_cycles['R3']),
        'index_length': read_length(*read_cycles['Ri']),
    }

    print("""\
[read_format]
total-cycles = {cv[total_cycles]}
fivep-start = {read_cycles[R5][0]}
fivep-length = {cv[fivep_length]}
index-start = {read_cycles[Ri][0]}
index-length = {cv[index_length]}
threep-start = {read_cycles[R3][0]}
threep-length = {cv[threep_length]}
""".format(cv=computed_values, read_cycles=read_cycles), file=outf)


def generate_options_section(outf):
    print("""\
[options]
keep-no-delimiter = {of[keep_no_delimiter]:d}
keep-low-quality-balancer = {of[keep_low_quality_balancer]:d}
threads = {threads}
read-buffer-size = {pf[maximum_buffer_size]}
""".format(of=params.conf['output_filtering'], pf=params.conf['performance'],
           threads=threads), file=outf)


def generate_output_section(outf):
    print("""\
[output]
seqqual = scratch/seqqual/{{name}}_{tile}.txt.gz
taginfo = scratch/taginfo/{{name}}_{tile}.txt.gz
stats = scratch/stats/signal-proc-{tile}.csv
length-dists = scratch/stats/length-dist-{tile}.csv
signal-dump = scratch/sigdumps/signaldump-{{name}}-{tile}.dmp.gz
""".format(tile=wildcards.tile), file=outf)

    for subdir in 'scratch/seqqual scratch/taginfo scratch/stats scratch/sigdumps'.split():
        if not os.path.isdir(subdir):
            os.makedirs(subdir)


def generate_alternative_calls_section(outf):
    print("[alternative_calls]", file=outf)
    for readname, program in params.conf['third_party_basecaller'].items():
        if program:
            startcycle = params.conf['read_cycles'][readname][0]
            seqpath = 'scratch/{prog}calls/{name}_{tile}.fastq.gz'.format(
                            prog=program.lower(), name=readname, tile=wildcards.tile)
            print('{} = {}'.format(startcycle, seqpath), file=outf)

    print('', file=outf)


def generate_control_section(outf):
    print("""\
[control]
phix-match-name = PhiX
phix-match-start = {conf[phix_match_start]}
phix-match-length = {conf[phix_match_length]}
""".format(conf=params.conf['control']), file=outf)


def generate_balancer_section(outf):
    print("""\
[balancer]
start = {conf[start]}
length = {conf[length]}
minimum-occurrence = {conf[minimum_occurrence]}
minimum-quality = {conf[minimum_quality]}
minimum-qcpass-percent = {conf[minimum_quality_percent]}
num-positive-samples = {conf[positive_signal_samples]}
num-negative-samples = {conf[negative_signal_samples]}
""".format(conf=params.conf['balancer']), file=outf)


def generate_polyA_finder_section(outf):
    print("""\
[polyA_finder]
weight-T = {conf[weights][T]}
weight-A = {conf[weights][A]}
weight-C = {conf[weights][C]}
weight-G = {conf[weights][G]}
weight-N = {conf[weights][N]}
minimum-polya-length = {conf[minimum_polya_length]}
maximum-modifications = {conf[maximum_modifications]}
signal-analysis-trigger = {conf[signal_analysis_trigger]}
naive-count-trigger = {conf[naive_count_trigger]}
""".format(conf=params.conf['polyA_finder']), file=outf)


def generate_polyA_ruler_section(outf):
    print("""\
[polyA_ruler]
dark-cycles-threshold = {conf[dark_cycles_threshold]}
maximum-dark-cycles = {conf[maximum_dark_cycles]}
polya-score-threshold = {conf[polya_score_threshold]}
downhill-extension-weight = {conf[downhill_extension_weight]}
t-intensity-k = {conf[t_intensity_k]}
t-intensity-center = {conf[t_intensity_center]}
""".format(conf=params.conf['polyA_ruler']), file=outf)


def generate_experimental_samples_section(outf):
    for samplename in params.exp_samples + params.spikein_samples:
        print("""\
[sample:{name}]
index = {sf[index]}
maximum-index-mismatch = {sf[maximum-index-mismatch]}
delimiter-seq = {sf[delimiter-seq]}
delimiter-start = {sf[delimiter-start]}
maximum-delimiter-mismatch = {sf[maximum-delimiter-mismatch]}\
""".format(name=samplename, sf={
            'index':
                params.conf['experimental_samples'][samplename]
                if samplename in params.conf['experimental_samples']
                else params.conf['spikein_samples'][samplename],
            'maximum-index-mismatch': params.conf['maximum_mismatches']['index'][samplename],
            'delimiter-seq': params.conf['delimiter'][samplename][1],
            'delimiter-start': params.conf['delimiter'][samplename][0],
            'maximum-delimiter-mismatch': params.conf['maximum_mismatches'
                                                      ]['delimiter'][samplename],
        }), file=outf)

        if samplename in params.conf['dupcheck_regions']:
            for i, (start, end) in enumerate(params.conf['dupcheck_regions'][samplename]):
                print("""\
umi-start:{num} = {start}
umi-length:{num} = {length}\
""".format(num=i + 1, start=start, length=end-start+1), file=outf)

        if samplename in params.conf['fingerprint_fixed_sequence']:
            print("""\
fingerprint-seq = {seq}
fingerprint-start = {start}
maximum-fingerprint-mismatch = {mismatches}\
""".format(seq=params.conf['fingerprint_fixed_sequence'][samplename],
           start=params.conf['read_cycles']['R3'][0], # TODO: make this configurable?
           mismatches=params.conf['maximum_mismatches']['fingerprint'][samplename]),
           file=outf)

        if samplename in params.conf['spikein_trimming_length']:
            print("limit-threep-processing =",
                  params.conf['spikein_trimming_length'][samplename], file=outf)

        if samplename in params.conf['debug']['dump_signals']:
            print("dump-processed-signals =",
                  int(params.conf['debug']['dump_signals'][samplename]), file=outf)

        print('', file=outf)


if is_snakemake_child:
    import sys

    outf = open(output.sigproc_conf, 'w')

    generate_source_section(outf)
    generate_read_format_section(outf)
    generate_options_section(outf)
    generate_output_section(outf)
    generate_alternative_calls_section(outf)
    generate_control_section(outf)
    generate_balancer_section(outf)
    generate_polyA_finder_section(outf)
    generate_polyA_ruler_section(outf)
    generate_experimental_samples_section(outf)

