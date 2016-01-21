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

from tailseeker.fileutils import LineParser, open_gzip_buffered
from base64 import b64decode
import numpy as np
import sys

if sys.version_info[0] >= 3:
    from tailseqext import decode_intensity
else:
    from tailseqext2 import decode_intensity

readid_key = lambda x: (x.tile, x.cluster)
decode_ascii = lambda x: x.decode('ascii')


def unpack_readid(s):
    tile, cluster = s.split(b':')
    return tile, int(cluster)

def pack_readid(tile, cluster):
    return '{}:{:08d}'.format(tile, cluster)

parse_readid_list = LineParser([
    ('tile', decode_ascii),
    ('cluster', int),
], linefeed=b'\r\n', separator=b':')

parse_sqi_lite = LineParser([
    ('tile', decode_ascii),
    ('cluster', int),
    ('istart', int),
    ('seq', decode_ascii),
    ('qual', None),
    ('intensity', None),
], linefeed=b'\r\n')

def decode_phred_quality(s, scale=33):
    return np.fromstring(s, np.uint8) - scale

parse_sqi = LineParser([
    ('tile', decode_ascii),
    ('cluster', int),
    ('istart', int), # 0-based coordinate of insert start position
    ('seq', decode_ascii),
    ('qual', decode_phred_quality),
    ('intensity', decode_intensity),
], linefeed=b'\r\n')

def decode_pascore(s):
    return np.fromstring(b64decode(s), np.float32)

parse_pascore = LineParser([
    ('tile', decode_ascii),
    ('cluster', int),
    ('endmod_len', int),
    ('seqbased_polya_len', int),
    ('pascore', decode_pascore),
], linefeed=b'\r\n')

parse_polya_calls = LineParser([
    ('tile', decode_ascii),
    ('cluster', int),
    ('start_pos', int), # 0-based, relative to the immediate next base to the delimiter
    ('polya_len', int),
    ('seqbased_len', int),
    ('hmmbased_len', int),
    ('naive_len', int),
], linefeed=b'\r\n')

