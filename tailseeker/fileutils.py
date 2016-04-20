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

from itertools import groupby, chain
from struct import unpack
import subprocess as sp
import random
import shutil
import time
import tempfile
import os
import sys
import glob
import gzip
import io


__all__ = [
    'LineParser', 'ParsedLine', 'ParsedLineComment', 'TemporaryDirectory',
    'ParallelMatchingFilter', 'ParallelMatchingReader',
    'open_gzip_pipe', 'open_gzip_buffered', 'MultiJoinIterator',
    'open_bgzip_writer', 'merge_bgzf_files',
]

BASH_CMD = os.environ.get('BASH_CMD', '/bin/bash')
BGZIP_CMD = os.environ.get('BGZIP_CMD', 'bgzip')


BGZF_EOF_BLOCK = (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03"
                  b"\x00\x00\x00\x00\x00\x00\x00\x00\x00")

class LineParser(object):

    def __init__(self, fields_spec, separator=b'\t', linefeed=b'\n', comment=None, listtrailer=False):
        self.field2index = dict((field[0], i) for i, field in enumerate(fields_spec))
        self.fields_spec = fields_spec
        self.separator = separator
        self.linefeed = linefeed
        self.comment = comment
        self.listtrailer = listtrailer

    def parse(self, line):
        if self.comment is not None and line.startswith(self.comment):
            return ParsedLineComment(line.rstrip(self.linefeed), line)

        fields = line.rstrip(self.linefeed).split(self.separator)
        if self.listtrailer:
            list_start_index = len(self.fields_spec) - 1
            fields = fields[:list_start_index] + [fields[list_start_index:]]

        return ParsedLine(self.field2index, [
            (data if adapter is None else (adapter(data) if data is not None else None))
            for data, (name, adapter) in zip(fields, self.fields_spec)], line)

    def iter_parse(self, it):
        for line in it:
            yield self.parse(line)

    def __call__(self, it):
        if isinstance(it, str):
            it = open_gzip_buffered(it)

        return self.iter_parse(it)

    def as_table(self, inputfile, **kwds):
        import pandas
        fieldnames = [fn for fn, _ in self.fields_spec]
        if isinstance(inputfile, str):
            inputfile = open_gzip_buffered(inputfile)
        return pandas.read_table(inputfile, names=fieldnames, **kwds)

    # for compatibility, May 29, 2014.
    to_table = as_table

class ParsedLineComment(object):

    def __init__(self, data, line):
        self.data = data
        self.line = line

    if sys.version_info[0] >= 3:
        def __str__(self):
            return self.line.decode()
    else:
        def __str__(self):
            return self.line


class ParsedLine(object):

    def __init__(self, field2index, data, line):
        self.field2index = field2index
        self.data = data
        self.line = line

    def __getitem__(self, index):
        return self.data[index]

    def __len__(self):
        return len(self.data)

    def __getattr__(self, name):
        if name.startswith('_'):
            return object.__getattr__(self, name)
        else:
            return self.data[self.field2index[name]]

    def __getitem__(self, name):
        return self.data.__getitem__(name)

    def __repr__(self):
        return '<ParsedLine %s>' % ' '.join(
            '%s=%s' % (name, repr(self.data[idx]))
            for name, idx in sorted(self.field2index.items(), key=lambda x: x[1])
            if idx < len(self.data))

    if sys.version_info[0] >= 3:
        def __str__(self):
            return self.line.decode()
    else:
        def __str__(self):
            return self.line


class ParallelMatchingReader(object):

    def __init__(self, src1, src2, src1key, src2key=None, separator='\t', linefeed='\n'):
        if src2key is None:
            src2key = src1key

        if callable(src1key) and callable(src2key):
            self.src1 = groupby(src1, key=src1key)
            self.src2 = groupby(src2, key=src2key)
        else:
            self.src1 = groupby(self.splititer(src1, separator, linefeed),
                                key=lambda x: x[src1key])
            self.src2 = groupby(self.splititer(src2, separator, linefeed),
                                key=lambda x: x[src2key])
            self.src1key = src1key
            self.src2key = src2key

    @staticmethod
    def splititer(src, sep, lf):
        for line in src:
            yield line.rstrip(lf).split(sep)

    def __iter__(self):
        src1ontap = src2ontap = None

        while True:
            if src1ontap is None:
                src1ontap = next(self.src1, None)
            if src2ontap is None:
                src2ontap = next(self.src2, None)

            if src1ontap is None:
                if src2ontap is None:
                    # both inputs are null
                    break

                yield (None, src2ontap[1])
                for _, v in self.src2:
                    yield (None, v)

                break
            elif src2ontap is None:
                yield (src1ontap[1], None)
                for _, v in self.src1:
                    yield (v, None)

                break

            # inputs from both sources are available
            if src1ontap[0] < src2ontap[0]:
                yield (src1ontap[1], None)
                src1ontap = None
            elif src1ontap[0] > src2ontap[0]:
                yield (None, src2ontap[1])
                src2ontap = None
            else:
                yield (src1ontap[1], src2ontap[1])
                src1ontap = src2ontap = None


class ParallelMatchingFilter(ParallelMatchingReader):

    def __iter__(self):
        baseit = super(ParallelMatchingFilter, self).__iter__()
        for matcha, matchb in baseit:
            if matcha is not None and matchb is not None:
                for item in matcha:
                    yield item


class _IteratorFinished(object):

    def false_funcs(self, other):
        return False
    __lt__ = __le__ = __eq__ = false_funcs

    def true_funcs(self, other):
        return True
    __gt__ = __ge__ = __ne__ = true_funcs

    def __nonzero__(self):
        return True

_IteratorFinished = _IteratorFinished()


class MultiJoinIterator(object):

    def __init__(self, sources, keyfuncs):
        if callable(keyfuncs):
            keyfuncs = [keyfuncs] * len(sources)

        self.groupsrcs = [groupby(src, keyfun)
                          for src, keyfun in zip(sources, keyfuncs)]

    def __iter__(self):
        nextitems = [None] * len(self.groupsrcs)

        while True:
            nextitems = [
                item if item is not None else next(self.groupsrcs[srcno], _IteratorFinished)
                for srcno, item in enumerate(nextitems)
            ]

            if all(item is _IteratorFinished for item in nextitems):
                break

            minkey = min(item[0] for item in nextitems if item is not _IteratorFinished)
            output = [minkey]
            for itemno, item in enumerate(nextitems):
                if item is _IteratorFinished or item[0] != minkey:
                    output.append(iter(''))
                else:
                    output.append(item[1])
                    nextitems[itemno] = None

            yield output


class TemporaryDirectory(object):
    def __init__(self, dir=None, asobj=False, automerge=False):
        self.dir = os.environ.get('TAILSEQ_SCRATCH_DIR', '.') if dir is None else dir
        self.path = None
        self.asobj = asobj
        self.output_counter = 0
        self.automerge = automerge

    def __enter__(self):
        self.path = tempfile.mkdtemp(dir=self.dir)
        return self if self.asobj else self.path

    def __exit__(self, type, value, traceback):
        if self.path is not None:
            if self.automerge:
                self.merge_into_stdout()
            shutil.rmtree(self.path)

    def __str__(self):
        return self.path

    def all_files(self):
        return sorted(glob.glob(os.path.join(self.path, '*')))

    def __iter__(self):
        return chain(*map(open, self.all_files()))

    def merge_into_file(self, outfile):
        infiles = self.all_files()
        if infiles:
            sp.check_call(['cat'] + infiles, stdout=outfile)

    def merge_into_stdout(self):
        return self.merge_into_file(sys.stdout)

    def merge_bgzf(self, outfile, prefix=None):
        files = self.all_files()
        if prefix is not None:
            files = [f for f in files
                     if os.path.basename(f).startswith(prefix)]
        merge_bgzf_files(outfile, files)

    def next_output_file(self):
        newpath = os.path.join(self.path, format(self.output_counter, '08x'))
        self.output_counter += 1
        return newpath


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


def merge_bgzf_files(outputfile, inputfiles):
    output = open(outputfile, 'wb')

    for inpfile in inputfiles:
        try:
            clone_bgzf_blocks(output, open(inpfile, 'rb'))
        except StopIteration:
            pass

    output.write(BGZF_EOF_BLOCK)


def open_gzip_pipe(filename):
    return sp.Popen(['zcat', filename], stdout=sp.PIPE).stdout


def open_gzip_buffered(filename, mode='rb'):
    if open(filename, 'rb').read(2) == b'\x1f\x8b':
        if 't' in mode:
            return io.TextIOWrapper(io.BufferedReader(gzip.open(filename)))
        else:
            return io.BufferedReader(gzip.open(filename))
    else:
        return open(filename, mode)


def open_bgzip_writer(filename, mode='b'):
    subproc = sp.Popen('"{bgzip}" -c /dev/stdin > "{output}"'.format(
                            bgzip=BGZIP_CMD, output=filename), shell=True, stdin=sp.PIPE)
    return (
        io.TextIOWrapper(subproc.stdin) if 't' in mode else subproc.stdin,
        subproc,
    )


if __name__ == '__main__':
    for a, b in ParallelMatchingReader(open('t2ctemp.bed'), open('new-t2c-genome.bed'), 3, 3):
        print(list(a or ''), list(b or ''))

