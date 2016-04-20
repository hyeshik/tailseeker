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
