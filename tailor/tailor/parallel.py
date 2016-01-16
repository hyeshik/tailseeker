#!/usr/bin/env python3
#
# Copyright (c) 2013 Institute for Basic Science
# Copyright (c) 2013 Hyeshik Chang
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
import random
import os

__all__ = ['open_tabix_parallel']

if 'TABIX_CMD' in os.environ:
    TABIX_CMD = os.environ['TABIX_CMD']
else:
    TABIX_CMD = 'tabix'


class TabixOpener(object): # must be pickleable so that can be passed through futures

    def __init__(self, bgzfile, interval):
        self.bgzfile = bgzfile
        self.interval = interval

    def __call__(self):
        proc = sp.Popen([TABIX_CMD, self.bgzfile, self.interval], stdout=sp.PIPE)
        return proc.stdout

    def __len__(self):
        return int(sp.check_output('"{}" "{}" "{}" | wc -l'.format(
                    TABIX_CMD, self.bgzfile, self.interval), shell=True).strip())

    def random_sample(self, num):
        # XXX: this can be improved using reservoir sampling (Vitter, 1995)
        total_records = len(self)
        sampled_rec_no = set(random.sample(xrange(total_records), num))

        for recno, line in enumerate(self()):
            if recno in sampled_rec_no:
                yield line


def open_tabix_parallel(bgzfile, startpos=0, endpos=(2**31-1), named=False):
    regions = sp.check_output([TABIX_CMD, '-l', bgzfile]).decode('ascii').split()

    if named:
        return dict((reg, TabixOpener(bgzfile, '{}:{}-{}'.format(reg, startpos, endpos)))
                    for reg in regions)
    else:
        return [TabixOpener(bgzfile, '{}:{}-{}'.format(reg, startpos, endpos))
                for reg in regions]


if __name__ == '__main__':
    import pickle
    print(pickle.dumps(open_tabix_parallel('sequences/third_118.sqi.bgz')))

