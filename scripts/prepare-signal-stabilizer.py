#!/usr/bin/env python3
#
# Copyright (c) 2014-5 Institute for Basic Science
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

import pickle
import time
import numpy as np
from tailor.signalproc import TAILseqSignalProcessor

def parse_interval(s):
    begin, end = map(int, s.split(':'))
    return slice(begin - 1, end)


def inspect_signals(options, openinput):
    from tailor.parsers import parse_sqi

    highinterval = parse_interval(options.proberange)
    highrange_width = highinterval.stop - highinterval.start

    HIGHPROBE_PERCENTILE = 10
    HIGHPROBE_MINIMUM_VALUE = 1.0

    readinterval = parse_interval(options.readrange)

    highinspection_interval = parse_interval(options.highinspection)
    highinspection_basis_interval = parse_interval(options.highinspection_basis)

    tilenumber = None
    signalproc = TAILseqSignalProcessor(options.cyclemtx, options.colormtx, options.readno,
                                        readinterval.start, readinterval.stop,
                                        options.spotnormlen, options.lowsignalmask)

    highprobed_scale_region = []
    highprobed_scale_basis = []

    REQUIRED_TERMINAL_T_STRETCH = 'T' * 3

    for spot in parse_sqi(openinput()):
        if tilenumber is None:
            tilenumber = spot.tile

        if (REQUIRED_TERMINAL_T_STRETCH and
                not spot.seq[readinterval][highinspection_interval].startswith(REQUIRED_TERMINAL_T_STRETCH)):
            continue

        t_stand_out = signalproc.calc_t_stand_out(tilenumber, spot.seq, spot.intensity)
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
    return {
        tilenumber:
        refbasis / np.median(np.array(highprobed_scale_region), axis=0)
    }

def run(options):
    from tailor.parallel import open_tabix_parallel, TabixOpener
    from concurrent import futures

    infile = options.infile[0]
    #print inspect_signals(options, open_tabix_parallel(infile)[0], options.output)
    scale_factors = {}
    configuration = {
        'ctime': time.time(),
        'readno': options.readno,
        'color_matrix': options.colormtx,
        'cycle_matrix': options.cyclemtx,
        'high_probe_pAref_range': parse_interval(options.proberange),
        'high_probe_scale_range': parse_interval(options.highinspection),
        'high_probe_basis_range': parse_interval(options.highinspection_basis),
        'read_range': parse_interval(options.readrange),
        'spot_norm_length': options.spotnormlen,
        'low_signal_mask': options.lowsignalmask,
    }

    if options.parallel <= 1:
        for jobno, opener in enumerate(open_tabix_parallel(infile)):
            scale_factors.update(inspect_signals(options, opener))

        pickle.dump((configuration, scale_factors), open(options.output, 'wb'))

        return

    with futures.ProcessPoolExecutor(options.parallel) as executor:
        jobs = []

        for jobno, opener in enumerate(open_tabix_parallel(infile)):
            job = executor.submit(inspect_signals, options, opener)
            jobs.append(job)

        for j in jobs:
            scale_factors.update(j.result())

        pickle.dump((configuration, scale_factors), open(options.output, 'wb'))


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Inspect T-stand-out signals from a sample'
                                                 ' and produce scaling factors')
    parser.add_argument(dest='infile', metavar='FILE', type=str, nargs=1,
                        help='Path to a sqi file')
    parser.add_argument('--parallel', dest='parallel', metavar='NUMBER', type=int, default=1,
                        help='Number of parallel processes')
    parser.add_argument('--output', dest='output', metavar='FILE', type=str, required=True,
                        help='Output file for collected statistics')
    parser.add_argument('--read', dest='readno', metavar='NUMBER', type=str, required=True,
                        help='Read number to process (3 in the standard paired-end setup).')
    parser.add_argument('--color-matrix', dest='colormtx', metavar='FILE', type=str,
                        required=True,
                        help='Path to a pickle file that contains color matrices.')
    parser.add_argument('--cycle-scaling', dest='cyclemtx', metavar='FILE', type=str,
                        required=True,
                        help='Path to a pickle file that contains per-scale normalization '
                             'matrices.')
    parser.add_argument('--high-probe-range', dest='proberange', metavar='START:END', type=str,
                        default='70:120',
                        help='Cycle number range in inserts for long poly(A) tail detection.')
    parser.add_argument('--high-probe-scale-inspection', dest='highinspection',
                        metavar='START:END', type=str,
                        default='21:70',
                        help='Cycle number range in inserts for scale factor inspection from '
                             'long poly(A) tails.')
    parser.add_argument('--high-probe-scale-basis', dest='highinspection_basis',
                        metavar='START:END', type=str, default='71:75',
                        help='Cycle number reange in inserts for basis of scale factors')
    parser.add_argument('--read-range', dest='readrange', metavar='start:end', type=str,
                        required=True,
                        help='Cycle number range of the read in 1-based inclusive coordinate.')
    parser.add_argument('--spot-norm-length', dest='spotnormlen', metavar='LEN', type=int,
                        default=20,
                        help='Use the first N cycles for detection of spot brightness')
    parser.add_argument('--low-signal-mask', dest='lowsignalmask', metavar='LIMIT',
                        type=float, default=0.3,
                        help='Mask low signals that have less sum of signals than this')
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    options = parse_arguments()
    run(options)

