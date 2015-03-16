#!/usr/bin/env python3
#
# Copyright (c) 2015 Hyeshik Chang
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

TILE_LIST = {
    'HiSeq-RapidV1': """
        1101 1102 1103 1104 1105 1106 1107 1108
        1109 1110 1111 1112 1113 1114 1115 1116
        1201 1202 1203 1204 1205 1206 1207 1208
        1209 1210 1211 1212 1213 1214 1215 1216
        2101 2102 2103 2104 2105 2106 2107 2108
        2109 2110 2111 2112 2113 2114 2115 2116
        2201 2202 2203 2204 2205 2206 2207 2208
        2209 2210 2211 2212 2213 2214 2215 2216""".split(),
    'MiSeq-V2': """
        1101 1102 1103 1104 1105 1106 1107 1108 1109 1110 1111 1112 1113 1114
        2101 2102 2103 2104 2105 2106 2107 2108 2109 2110 2111 2112 2113 2114""".split(),
}

SIGNAL_SCALES = {
    'HiSeq-RapidV1': 2,
    'MiSeq-V2': 0,
}


def get_tiles(conf):
    tilemaps = {}

    for source in conf['sources']:
        try:
            tiles = TILE_LIST[source['type']]
        except KeyError:
            raise ValueError("{} is in the list of known sequencing chemistry.".format(
                                source['type']))

        for tile in tiles:
            tileid = source['id'] + tile
            tilemaps[tileid] = {
                'datadir': source['dir'],
                'lane': source['lane'],
                'type': source['type'],
                'tile': tile,
                'id': tileid,
                'laneid': source['id'],
            }

    return tilemaps

def get_signalscale(chemistry_type):
    return SIGNAL_SCALES[chemistry_type]

