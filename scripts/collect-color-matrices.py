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

import numpy as np
import os
import re
import pickle
import base64

def load_colormatrix(filename):
    elements = list(map(float, open(filename).read().split()))
    return np.array(elements).reshape(4, 4)

def load_decrosstalk_matrices(tilemapping):
    matrices = {}

    for tile, path in tilemapping.items():
        tile = str(tile)
        path = str(path)

        allfiles = os.listdir(os.path.dirname(path))
        fnpattern = re.compile(os.path.basename(path).replace('READNO', '(\d+)'))
        for match in filter(None, (fnpattern.search(fn) for fn in allfiles)):
            readno = match.groups()[0]
            fullfn = os.path.join(os.path.dirname(path), match.group())
            matrices[tile, readno] = load_colormatrix(fullfn)

    return matrices

def run(options):
    tilemapping = pickle.loads(base64.b64decode(options.tilemapping))
    matrices = load_decrosstalk_matrices(tilemapping)
    pickle.dump(matrices, open(options.output, 'wb'))

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Collect color matrices from the original '
                                                 'run output')
    parser.add_argument('--tile-mapping', dest='tilemapping', metavar='STR', type=str,
                        required=True,
                        help='base64-encoded pickle dump of a dictionary')
    parser.add_argument('--output', dest='output', metavar='FILE', type=str, required=True,
                        help='Output file for collected matrices')
    options = parser.parse_args()

    return options

if __name__ == '__main__':
    options = parse_arguments()
    run(options)

