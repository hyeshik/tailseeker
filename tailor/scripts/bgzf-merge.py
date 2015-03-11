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

import subprocess as sp
from struct import unpack


BGZF_EOF_BLOCK = (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03"
                  b"\x00\x00\x00\x00\x00\x00\x00\x00\x00")


def clone_bgzf_blocks(output, input):
    while True:
        header = input.read(12)
        if not header:
            raise StopIteration

        h_id, h_cm, h_flg, _, _, _, h_extra_len = unpack(b'<HBBLBBH', header)
        if (h_id, h_cm, h_flg) != (0x8b1f, 8, 4): # the fixed magic of BGZF
            raise ValueError("Not a regular BGZF file.")

        extra_block = input.read(h_extra_len)
        subblk_start = 0
        block_size = None
        while subblk_start < len(extra_block):
            subblk_id = extra_block[subblk_start:subblk_start+2]
            subblk_len = unpack(b'<H', extra_block[subblk_start+2:subblk_start+4])[0]
            if subblk_id == b'BC': # BGZF ext field id
                subblk_data = extra_block[subblk_start+4:subblk_start+4+subblk_len]
                block_size = unpack(b'<H', subblk_data)[0]
            subblk_start += subblk_len + 4

        if block_size is None:
            raise ValueError("A mandatory BGZF extra field is not found.")

        rest_of_block = input.read(block_size - h_extra_len - 11)
        input_size = rest_of_block[-4:]

        if input_size != b'\x00\x00\x00\x00': # not an EOF block
            output.write(header)
            output.write(extra_block)
            output.write(rest_of_block)


def process(options):
    output = open(options.output, 'wb')

    for inpfile in options.inputs:
        try:
            clone_bgzf_blocks(output, open(inpfile, 'rb'))
        except StopIteration:
            pass

    output.write(BGZF_EOF_BLOCK)


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
    options = parse_arguments()
    process(options)

