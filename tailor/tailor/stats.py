#!/usr/bin/env python3
#
# Copyright (c) 2013-2015 Institute for Basic Science
# Copyright (c) 2012-2013 Hyeshik Chang
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
    'similarity_sort', 'sample_iterable', 'ReservoirSampler', 'smooth',
    'savitzky_golay',
]

from Bio.Cluster import cluster
from scipy.spatial.distance import pdist, squareform
import random
import numpy as np

def similarity_sort(data, dist='correlation'):
    #tree = cluster.treecluster(data, dist='c')
    if isinstance(dist, str):
        distmatrix = squareform(pdist(data, dist))
    else:
        distmatrix = dist
    tree = cluster.treecluster(distancematrix=distmatrix)
    pushed = {}
    def resolve(idx):
        if idx < 0:
            getting = pushed[idx]
            del pushed[idx]
            return getting
        else:
            return [idx]

    for i, node in enumerate(tree):
        clsi = -i-1
        pushed[clsi] = resolve(node.left) + resolve(node.right)

    return next(iter(pushed.values()))


class ReservoirSampler(object):

    def __init__(self, num):
        self.current = [None] * num
        self.num = num
        self.serial = 0

    def update(self, value, randrange=random.randrange):
        if self.serial < self.num:
            self.current[self.serial] = value
        else:
            s = randrange(self.serial)
            if s < self.num:
                self.current[s] = value

        self.serial += 1

    def __len__(self):
        return min(self.num, self.serial)

    def __getitem__(self, index):
        return self.current[index]

    def get(self):
        return self.current[:self.num]



def sample_iterable(iterable, num, randrange=random.randrange):
    # faster implementation of ReservoirSampler

    result = [item for i, item in zip(range(num), iterable)]
    if len(result) < num:
        return result

    for i, item in enumerate(iterable, num):
        s = randrange(i)
        if s < num:
            result[s] = item

    return result


def group_consecutive(iterable, width=1):
    it = iter(iterable)
    first = next(it, None)
    if first is None:
        raise StopIteration

    stacked = [first]

    try:
        while True:
            nextvalue = next(it)
            if nextvalue - stacked[-1] <= width:
                stacked.append(nextvalue)
            else:
                yield stacked
                stacked = [nextvalue]

    except StopIteration:
        if stacked:
            yield stacked


# from http://wiki.scipy.org/Cookbook/SignalSmooth
def smooth(x, window_len=11, window='hanning'):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(),s,mode='same')

    return y[window_len:-window_len+1]

# from http://wiki.scipy.org/Cookbook/SavitzkyGolay
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

