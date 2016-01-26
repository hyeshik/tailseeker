#!/usr/bin/env python3
#
# Copyright (c) 2014-6 Institute for Basic Science
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
from tailseeker.signalproc import TAILseqSignalProcessor
from tailseeker.sequencers import tile_range_cond
import numpy as np
import tables
import pickle
import time


def inspect_signals(options, tileid):
    inputtablef = tables.open_file(options['infile'], 'r')
    nodegrp = '/sqi/{o[sample_name]}'.format(o=options)
    select_tile = tile_range_cond(tileid)
    seqquals = inputtablef.get_node(nodegrp + '/seqqual')
    intensities = inputtablef.get_node(nodegrp + '/intensities')

    highinterval = slice(*options['probe_range'])
    highrange_width = highinterval.stop - highinterval.start

    HIGHPROBE_PERCENTILE = 10
    HIGHPROBE_MINIMUM_VALUE = 1.0

    readinterval = slice(*options['read_range'])

    highinspection_interval = slice(*options['high_inspection'])
    highinspection_basis_interval = slice(*options['high_basis'])

    signalproc = TAILseqSignalProcessor(options['global_scaling'],
                                        options['color_matrix'],
                                        options['readno'],
                                        readinterval.start, readinterval.stop,
                                        options['spot_norm_length'],
                                        options['low_signal_mask'])

    highprobed_scale_region = []
    highprobed_scale_basis = []

    REQUIRED_TERMINAL_T_STRETCH = b'TTT'

    for sqrow in seqquals.where(select_tile):
        introw = intensities[sqrow.nrow]
        seq = sqrow['seq']

        if (REQUIRED_TERMINAL_T_STRETCH and
                not seq[readinterval][highinspection_interval].startswith(REQUIRED_TERMINAL_T_STRETCH)):
            continue

        t_stand_out = signalproc.calc_t_stand_out(tileid, seq, introw)
        if t_stand_out is None:
            continue

        #print np.round(t_stand_out, 2)
        highprobe_region = t_stand_out[highinterval]
        if (len(highprobe_region) >= highrange_width and
                    np.percentile(highprobe_region, HIGHPROBE_PERCENTILE)
                            > HIGHPROBE_MINIMUM_VALUE):
            highprobed_scale_region.append(t_stand_out[highinspection_interval])
            highprobed_scale_basis.append(t_stand_out[highinspection_basis_interval])

    refbasis = np.median(np.array(highprobed_scale_basis), axis=0).mean()

    return {tileid: refbasis / np.median(np.array(highprobed_scale_region), axis=0)}


from tailseeker.sequencers import decode_gcid, tile_range_cond
def get_tiles(seqquals):
    tileid_lo = decode_gcid(seqquals[0]['gcid'])[0]
    tileid_hi = decode_gcid(seqquals[-1]['gcid'])[0]
    return list(range(tileid_lo, tileid_hi + 1))


def write_out_signal_parameters(options, iter_results):
    with tables.open_file(options['output'], 'w') as outf:
        sigprocgrp = outf.create_group('/', 'signalproc',
                                       'Parameters related to signal processing')
        samplegrp = outf.create_group('/signalproc', 'persample',
                                      'Signal scale normalization parameters for each sample')
        thissamplegrp = outf.create_group(samplegrp, options['sample_name'],
                                          'Scaling normalization ranges')

        outf.create_array(thissamplegrp, 'conf', [pickle.dumps(options)])

        for tileid, scalingarray in iter_results:
            tilename = 'tile_{}'.format(tileid)
            outf.create_array(thissamplegrp, tilename, scalingarray)


def run(options):
    with tables.open_file(options['infile'], 'r') as tablef:
        seqquals = tablef.get_node('/sqi/{o[sample_name]}/seqqual'.format(o=options))
        tiles = get_tiles(seqquals)

    scale_factors = {}
    options['ctime'] = time.time()

    if options['parallel'] <= 1:
        for tileid in tiles:
            scale_factors.update(inspect_signals(options, tileid))
    else:
        from concurrent import futures

        with futures.ProcessPoolExecutor(options['parallel']) as executor:
            jobs = []

            for tileid in tiles:
                job = executor.submit(inspect_signals, options, tileid)
                jobs.append(job)

            for j in jobs:
                scale_factors.update(j.result())

    write_out_signal_parameters(options, scale_factors.items())


if is_snakemake_child:
    convert_1based = lambda x: (x[0] - 1, x[1])
    convert_r3_range = lambda x: (balancer_length + x[0] - 1, balancer_length + x[1])

    options = {
        'parallel': threads,
        'infile': input.signals,
        'output': output[0],
        'readno': read_cycles['R3'][2],
        'sample_name': wildcards.sample,
        'color_matrix': input.colormatrix,
        'global_scaling': input.globalscaling,
        'probe_range': convert_r3_range(normalizer_conf['polya_detection_range']),
        'high_inspection': convert_r3_range(normalizer_conf['target_range']),
        'high_basis': convert_r3_range(normalizer_conf['reference_range']),
        'read_range': convert_1based(read_cycles['R3'][:2]),
        'spot_norm_length': balancer_length,
        'low_signal_mask': normalizer_conf['low_signal_mask'],
    }

    run(options)
