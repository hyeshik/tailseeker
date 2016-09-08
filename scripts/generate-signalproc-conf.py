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
signal = scratch/signals/{{name}}_{tile}.sigpack
signal-dists = scratch/sigdists-r00/{{posneg}}_{tile}.sigdists
stats = scratch/stats/signal-proc-{tile}.csv
length-dists = scratch/stats/length-dist-{tile}.csv
""".format(tile=wildcards.tile), file=outf)

    for subdir in 'seqqual taginfo signals sigdists stats'.split():
        subdir = os.path.join('scratch', subdir)
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


def generate_polyA_seeder_section(outf):
    print("""\
[polyA_seeder]
seed-trigger-polya-length = {conf[seed_trigger_polya_length]}
negative-sample-polya-length = {conf[negative_sample_polya_length]}
max-cctr-scan-left-space = {conf[max_cctr_scan_left_space]}
max-cctr-scan-right-space = {conf[max_cctr_scan_right_space]}
required-cdf-contrast = {conf[required_cdf_contrast]}
polya-boundary-pos = {conf[polya_boundary_pos]}
polya-sampling-gap = {conf[polya_sampling_gap]}
dist-sampling-bins = {conf[dist_sampling_bins]}
fair-sampling-fingerprint-length = {conf[fair_sampling_fingerprint_length]}
fair-sampling-hash-space-size = {conf[fair_sampling_hash_space_size]}
fair-sampling-max-count = {conf[fair_sampling_max_count]}
""".format(conf=params.conf['polyA_seeder']), file=outf)


def generate_polyA_finder_section(outf):
    print("""\
[polyA_finder]
polyA-weight-T = {conf[polyA_weights][T]}
polyA-weight-A = {conf[polyA_weights][A]}
polyA-weight-C = {conf[polyA_weights][C]}
polyA-weight-G = {conf[polyA_weights][G]}
polyA-weight-N = {conf[polyA_weights][N]}
nonA-weight-T = {conf[nonA_weights][T]}
nonA-weight-A = {conf[nonA_weights][A]}
nonA-weight-C = {conf[nonA_weights][C]}
nonA-weight-G = {conf[nonA_weights][G]}
nonA-weight-N = {conf[nonA_weights][N]}
minimum-polya-length = {conf[minimum_polya_length]}
maximum-modifications = {conf[maximum_modifications]}
signal-analysis-trigger = {conf[signal_analysis_trigger]}
""".format(conf=params.conf['polyA_finder']), file=outf)


def generate_polyA_ruler_section(outf):
    print("""\
[polyA_ruler]
dark-cycles-threshold = {conf[dark_cycles_threshold]}
maximum-dark-cycles = {conf[maximum_dark_cycles]}
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
    generate_polyA_seeder_section(outf)
    generate_polyA_finder_section(outf)
    generate_polyA_ruler_section(outf)
    generate_experimental_samples_section(outf)

