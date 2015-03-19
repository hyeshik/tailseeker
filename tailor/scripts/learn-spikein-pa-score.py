#!/usr/bin/env python2
#
# Copyright (c) 2013 Institute for Basic Science
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

import ghmm
import numpy as np
import pickle
import random

SENTINEL_EMISSION = 1000000

DEFAULT_TRANSITION_MATRIX = np.array([
    #  A    B    M    N    S
    [ .94, .03, .01, .01, .01 ], # A 
    [ .0 , .5 , .4 , .08, .02 ], # B 
    [ .0 , .0 , .6 , .38, .02 ], # M 
    [ .0 , .0 , .0 , .95, .05 ], # N 
    [ .95, .01, .01, .03, .0  ], # S 

])

DEFAULT_EMISSION_GAUSSIANS_V1 = [ # (mu, sigma, weight)
    [(1.5, -1.0), (1.5, 1.5), (0.95, 0.05)], # A 
    [(1.5, -1.0), (1.5, 1.5), (0.75, 0.25)], # B 
    [(1.5, -1.0), (1.5, 1.5), (0.50, 0.50)], # M 
    [(1.5, -1.0), (1.5, 1.5), (0.25, 0.75)], # N 
    [(SENTINEL_EMISSION, SENTINEL_EMISSION), (0.5, 0.5), (0.5, 0.5)], # S 
]

DEFAULT_EMISSION_GAUSSIANS_V2 = [ # (mu, sigma, weight)
    [(1.5, 0.3), (.5, .4), (0.95, 0.05)], # A 
    [(1.5, 0.3), (.5, .4), (0.75, 0.25)], # B 
    [(1.5, 0.3), (.5, .4), (0.50, 0.50)], # M 
    [(1.5, 0.3), (.5, .4), (0.25, 0.75)], # N 
    [(SENTINEL_EMISSION, SENTINEL_EMISSION), (0.5, 0.5), (0.5, 0.5)], # S 
]

DEFAULT_START_PROBABILITIES = [
    0.0, # A 
    0.0, # B 
    0.0, # M 
    0.0, # N 
    1.0, # S
]

PARAMETER_PRESETS = {
    'v1': (DEFAULT_TRANSITION_MATRIX,
           DEFAULT_EMISSION_GAUSSIANS_V1,
           DEFAULT_START_PROBABILITIES),
    'v2': (DEFAULT_TRANSITION_MATRIX,
           DEFAULT_EMISSION_GAUSSIANS_V2,
           DEFAULT_START_PROBABILITIES),
}


F = ghmm.Float()


def learn(options):
    params = PARAMETER_PRESETS[options.preset]
    model = ghmm.HMMFromMatrices(F, ghmm.GaussianMixtureDistribution(F), *params)

    # load observations
    samplesets = []
    for inpname in options.input[0]:
        scorearr = np.load(inpname).clip(options.clipmin, options.clipmax)
        samplesets.extend([map(float, arr) for arr in scorearr])

    # randomize and flatten the observations
    random.shuffle(samplesets)

    trainingdata = [SENTINEL_EMISSION]
    for obs in samplesets:
        trainingdata.extend(obs)
        trainingdata.append(SENTINEL_EMISSION)

    trainingset = ghmm.EmissionSequence(F, trainingdata)

    # train the model
    model.baumWelch(trainingset, options.maxiter, 0.00001)

    return {
        'type': 'GMHMM',
        'score_version': (1 if options.preset == 'v1' else 2),
        'vrange': (options.clipmin, options.clipmax),
        'maxiter': options.maxiter,
        'parameters': model.asMatrices(),
    }


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Train a GMHMM from PA scores from '
                                                 'spike-in clusters')
    parser.add_argument(dest='input', metavar='FILE', action='append', nargs='+',
                        default=[], help='Sample PA score arrays to learn')
    parser.add_argument('--output', dest='output', metavar='FILE', required=True,
                        help='Output file for parameters of final model')
    parser.add_argument('--maximum-iterations', dest='maxiter', metavar='N', default=1000,
                        type=int, help='Maximum iteration count for Baum-Welch')
    parser.add_argument('--clip-minimum', dest='clipmin', metavar='MIN', default=-5.0,
                        type=float, help='Minimum value for input value clipping')
    parser.add_argument('--clip-maximum', dest='clipmax', metavar='MAX', default=5.0,
                        type=float, help='Maximum value for input value clipping')
    parser.add_argument('--preset', dest='preset', metavar='NAME', default='v1',
                        type=str, help='Initial parameters for the model')

    options = parser.parse_args()

    return options


if __name__ == '__main__':
    options = parse_arguments()

    model = learn(options)

    pickle.dump(model, open(options.output, 'wb'))

