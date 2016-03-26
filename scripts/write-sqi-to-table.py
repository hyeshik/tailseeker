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

from tailseeker.parsers import parse_sqi
from tailseeker.sequencers import encode_gcid
import numpy as np
import tables


def create_tables(h5, group, samplename, tileid, ncycles, nchannels):

    class SeqQual(tables.IsDescription):
        gcid = tables.Int64Col(pos=0)
        istart = tables.Int16Col(pos=1)
        seq = tables.StringCol(ncycles, pos=2)
        qual = tables.UInt8Col(shape=(ncycles,), pos=3)
        #intensity = tables.Int16Col(shape=(ncycles, nchannels), pos=5)

    subgroup = h5.create_group(group, samplename,
                    'SQI arrays for sample {} from tile {}'.format(samplename, tileid))

    seqqual = h5.create_table(subgroup, 'seqqual', SeqQual,
                              'Sequences and quality scores for sample {} from '
                              'tile {}'.format(samplename, tileid))
    intensities = h5.create_earray(subgroup, 'intensities', tables.Int16Atom(),
                                   (0, ncycles, nchannels),
                'Signal intensities for sample {} from tile {}'.format(samplename, tileid))

    return seqqual, intensities


def sqi2tbl(output_file, input_file, samplename, tileid, compression_level):
    filters = tables.Filters(complib='blosc:zlib', complevel=compression_level)

    with tables.open_file(output_file, mode='w',
                          title='Tailseeker seq-qual-intensity file',
                          filters=filters) as h5out:
        group = h5out.create_group('/', 'primary_sqi',
                                   'SQI arrays right after demultiplexing')
        seqqual = intensities = None

        for spot in parse_sqi(input_file):
            if seqqual is None: # the first row
                ncycles, nchannels = spot.intensity.shape # shape of intensity array
                seqqual, intensities = (
                    create_tables(h5out, group, samplename, tileid, ncycles, nchannels))

            gcid = encode_gcid(tileid, spot.cluster)

            seqqual.append([(gcid, spot.istart, spot.seq, spot.qual)])
            intensities.append([spot.intensity])

        if seqqual is not None:
            seqqual.cols.gcid.create_csindex()


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Convert a Sequence-Quality-Intensity file '
                                                 'input from stdin into a HDF5 storage.')
    parser.add_argument('--output', dest='output', metavar='FILE', type=str,
                        required=True, help='Path to the new file is located')
    parser.add_argument('--tileid', dest='tileid', metavar='NAME', type=int,
                        required=True, help='ID number of the tile currently processing')
    parser.add_argument('--samplename', dest='samplename', metavar='NAME', type=str,
                        required=True, help='Name of the sample currently processing')
    parser.add_argument('--compression-level', dest='comp_level', metavar='NUM', type=int,
                        default=5, help='Compression level for the output file')

    return parser.parse_args()


if __name__ == '__main__':
    import os

    options = parse_arguments()

    output_file = (
        options.output.replace('_SAMPLE_', options.samplename)
        if '_SAMPLE_' in options.output
        else options.output)

    sqi2tbl(output_file, os.fdopen(0, 'rb'),
            options.samplename, options.tileid,
            options.comp_level)

