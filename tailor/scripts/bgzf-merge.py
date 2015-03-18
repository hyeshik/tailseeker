#!/usr/bin/env python3
#
# Copyright (c) 2014 Institute for Basic Science
# Copyright (c) 2014 Hyeshik Chang
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


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Merges multiple bgzf files into one.')
    parser.add_argument('inputs', metavar='BGZF', type=str, nargs='+',
                        help='Input BGZF files.')
    parser.add_argument('--output', dest='output', metavar='OUTPUT', type=str,
                        help='Merged output BGZF file.', required=True)
    options = parser.parse_args()

    return options


if __name__ == '__main__':
    from tailor.fileutils import merge_bgzf_files

    options = parse_arguments()
    merge_bgzf_files(options.output, options.inputs)

