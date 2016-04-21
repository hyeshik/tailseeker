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

import numpy as np


PAFLAG_HAVE_3P_MODIFICATION         = 0x0001
PAFLAG_MEASURED_FROM_FLUORESCENCE   = 0x0002
PAFLAG_MEASURED_USING_NAIVE_RULER   = 0x0004
PAFLAG_DELIMITER_IS_SHIFTED         = 0x0008
PAFLAG_BARCODE_HAS_MISMATCHES       = 0x0010
PAFLAG_DARKCYCLE_EXISTS             = 0x0020
PAFLAG_DELIMITER_NOT_FOUND          = 0x0040
PAFLAG_DELIMITER_HAS_MISMATCH       = 0x0080
PAFLAG_BALANCER_CALL_QUALITY_BAD    = 0x0100
PAFLAG_BALANCER_BIASED              = 0x0200
PAFLAG_BALANCER_SIGNAL_BAD          = 0x0400
PAFLAG_DARKCYCLE_OVER_THRESHOLD     = 0x0800
PAFLAG_NO_POLYA_DETECTED            = 0x1000
PAFLAG_CHIMERIC_TEMPLATE            = 0x2000  # used by level 2 analysis
PAFLAG_LIKELY_HAVE_INTACT_END       = 0x4000  # used by level 2 analysis


seqid_format = '{r.tile}:{r.cluster:08d}:{r.pflags:04x}:{r.clones}:{r.polyA}:{r.mods}'

taginfo = {
    'compression': 'gzip',
    'names': ['tile', 'cluster', 'pflags', 'polyA', 'mods', 'clones'],
    'dtype': {'tile': str, 'cluster': np.uint32, 'pflags': np.uint32,
              'polyA': np.int32, 'mods': str, 'clones': np.int32},
    'keep_default_na': False,
}

refined_taginfo = {
    'compression': 'gzip',
    'names': ['tile', 'cluster', 'pflags', 'clones',
              'polyA', 'U', 'G', 'C', 'mods'],
    'dtype': {'tile': str, 'cluster': np.uint32, 'pflags': np.uint32,
              'clones': np.int32, 'polyA': np.int32, 'U': np.int32, 'G': np.int32,
              'C': np.int32, 'mods': str},
    'keep_default_na': False,
}
