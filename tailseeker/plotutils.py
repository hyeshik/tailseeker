#!/usr/bin/env python3
#
# Copyright (c) 2013 Institute for Basic Science
# Copyright (c) 2011-2013 Hyeshik Chang
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

__all__ = [
    'prepare_cumulative',
    'colormap_lch',
    'estimate_2d_density',
]

import numpy as np
from itertools import cycle
from colormath.color_objects import LCHabColor, sRGBColor
from colormath.color_conversions import convert_color

def colormap_lch(n, lum=(35, 65), chroma=75, start=0, end=300):
    if isinstance(lum, list) or isinstance(lum, tuple):
        lum = cycle(lum)
    else:
        lum = cycle([lum])

    rgbs = []
    for hue, lumn in zip(np.linspace(start, end, n), lum):
        rgb = convert_color(LCHabColor(lumn, chroma, hue), sRGBColor).get_value_tuple()
        rgbs.append('#{:02x}{:02x}{:02x}'.format(*map(int, np.array(rgb).clip(0, 1) * 255.0)))
    return rgbs

def prepare_cumulative(grp, width=100, reverse=False):
    grp = sorted(grp)
    xvalues = [grp[int(i / width * (len(grp)-1))] for i in range(width+1)]
    yvalues = np.arange(width+1) / width
    if reverse:
        yvalues = yvalues[::-1]
    return xvalues, yvalues

def estimate_2d_density(x, y, bw=None):
    from scipy.stats import gaussian_kde
    positions = np.vstack([x, y])
    kernel = gaussian_kde(positions, bw_method=bw)
    return kernel(positions)

