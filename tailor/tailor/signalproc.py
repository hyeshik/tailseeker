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
from numpy import linalg
from scipy.interpolate import interp1d

BASE2I = dict((b, i) for i, b in enumerate('ACGT'))
I2BASE = 'ACGT'
VERYSMALLNUMBER = 1e-300

# signal processor version 2
class TAILseqSignalProcessor:

    def __init__(self, normparamsfile, colormtxfile, readname, intv_start, intv_stop,
                       intensity_norm_length, lowsignalmask):
        self.intv_start = intv_start
        self.intv_stop = intv_stop
        self.read_interval = slice(intv_start, intv_stop)
        self.intensity_norm_length = intensity_norm_length

        scale_orig = pickle.load(open(normparamsfile, 'rb'))
        self.scaleparams = self.adopt_scale_params(scale_orig)
        self.colormtx = self.load_decrosstalk_matrices(colormtxfile, readname)

        self.lowsignalmask = lowsignalmask

    def load_decrosstalk_matrices(self, filename, readname):
        colormtx = pickle.load(open(filename, 'rb'))

        decrosstalk_matrices = {}
        for (vtile, readno), mtx in colormtx.items():
            if readno == readname:
                decrosstalk_matrices[vtile] = linalg.inv(mtx.T)

        return decrosstalk_matrices

    def adopt_scale_params(self, original_params):
        scale_params = {}

        for tile, params in original_params.items():
            structured_params = np.array([[params[cycle, base]
                                           for baseno, base in enumerate('ACGT')]
                                         for cycle in range(self.intv_start, self.intv_stop)])
            scale_params[tile] = structured_params.transpose()
            # shape = (2, 4, 251) in the standard setup

        return scale_params

    def calc_intensity_normalization_params(self, seq, intensity):
        seqarr = np.fromstring(seq, 'c')
        signal_params = []

        for basei, base in enumerate([b'A', b'C', b'G', b'T']):
            callpositions = np.where(seqarr == base)[0]
            if len(callpositions) < 1:
                return None

            posmean = intensity[callpositions, basei].mean()

            noncallpositions = np.where(seqarr != base)[0]
            if len(noncallpositions) < 1:
                return None
            negmean = intensity[noncallpositions, basei].mean()

            signal_params.append([negmean, posmean - negmean])

        return signal_params

    def calc_t_stand_out(self, tile, seq, intensity):
        # orthogonalize the signals
        adjusted_intensity = np.dot(intensity[self.read_interval], self.colormtx[tile])

        # calculate spot intensity normalization factors to fit into [0, 1] normally.
        intensity_norm_params = self.calc_intensity_normalization_params(
                seq[self.intv_start:self.intv_start+self.intensity_norm_length],
                adjusted_intensity[:self.intensity_norm_length])

        if intensity_norm_params is None: # not normalizable
            return None

        # (in-spot) normalize the intensity with the factors calculated above.
        normalized_intensity = np.array(
                [(adjusted_intensity[:, basei] - basalvalue) / dynrange
                for basei, (basalvalue, dynrange) in enumerate(intensity_norm_params)])

        # perform the per-cycle signal scaling
        percyclescaled_intensity = normalized_intensity.copy()
        for chan in range(4):
            percyclescaled_intensity[chan] -= self.scaleparams[tile][0, chan]
            percyclescaled_intensity[chan] /= self.scaleparams[tile][1, chan]

        # calculate T-(A+C+G) signal
        t_signal_stand_out = (percyclescaled_intensity[3]
                                    - percyclescaled_intensity[:3].sum(axis=0))

        # mask low signals (e.g. no reaction is proceeded further)
        darkcycles = np.where(percyclescaled_intensity.clip(0).sum(axis=0)
                                < self.lowsignalmask)[0]
        if len(darkcycles) > 0:
            t_signal_stand_out[darkcycles] = np.nan
            return self.fill_and_trim_nan(t_signal_stand_out)
        else:
            return t_signal_stand_out

    def fill_and_trim_nan(self, signals):
        nancnt = numcnt = 0
        isnan = np.isnan(signals)

        for i in range(len(signals) - 1, -1, -1):
            if isnan[i]:
                nancnt += 1
            else:
                numcnt += 1
                if numcnt * 5 > nancnt:
                    break

        for j in range(i, len(signals)):
            if isnan[j]:
                break

        signals = signals[:j]
        isnan = isnan[:j]

        if sum(isnan) == 0:
            return signals

        # impute erroneous spots using linear interpolation
        subindex = np.arange(len(signals))
        xvalues = subindex[~isnan]
        yvalues = signals[~isnan]
        if len(xvalues) < 3: # light spots are not enough.
            return None

        # out-of-bounds errors occur when first cycles are dark, thus fill with the first
        # available signal.
        fitfunc = interp1d(xvalues, yvalues, bounds_error=False, fill_value=yvalues[0],
                           kind='linear')

        signals[subindex[isnan]] = fitfunc(subindex[isnan])
        return signals

