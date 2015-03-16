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

from tailor.fileutils import LineParser, open_gzip_buffered
import numpy as np


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
