#!/usr/bin/env python3
#
# powersnake.py
#  - A set of common utility functions which are very useful in most Snakefiles.
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

__all__ = ['init_powersnake', 'external_script', 'init_powersnake',
           'load_snakemake_params', 'is_snakemake_child', 'suffix_filter',
           'tmpfile', 'notify']

from snakemake.shell import shell
import os
import tempfile

PARAMETER_PASSING_ENVVAR = 'SNAKEMAKE_PARAMS'
PUSHBULLET_ENVVAR = 'SNAKEMAKE_PUSHBULLET_APIKEY'


def is_snakemake_running():
    import inspect

    for outerframe in inspect.getouterframes(inspect.currentframe()):
        if outerframe[1].endswith('snakemake/__init__.py'):
            return True
    else:
        return False

is_snakemake_child = lambda: False

def init_powersnake():
    # pipefail is supported by bash only.
    shell.executable(os.popen('which bash').read().strip())
    shell.prefix('set -e; set -o pipefail; ')

def load_snakemake_params():
    global is_snakemake_child

    if PARAMETER_PASSING_ENVVAR not in os.environ:
        return

    is_snakemake_child = lambda: True

    import json, collections
    from snakemake.io import Namedlist
    import builtins

    options = json.load(open(os.environ[PARAMETER_PASSING_ENVVAR]))

    # Unset the parameter file var not to pass to children.
    del os.environ[PARAMETER_PASSING_ENVVAR]

    for varname, value in options.items():
        if not isinstance(value, int):
            value = Namedlist(fromdict=dict(value))
        setattr(builtins, varname, value)


def external_script(_command):
    import inspect, json, tempfile

    VARS_TO_PASS = 'input output threads wildcards params'.split()

    callerlocal = inspect.currentframe().f_back.f_locals
    callerglobal = inspect.currentframe().f_back.f_globals
    packed = {}
    for var in VARS_TO_PASS:
        if isinstance(callerlocal[var], int):
            packed[var] = callerlocal[var]
        else:
            numkeys = len(callerlocal[var].keys())
            if numkeys == 0: # when there is no input name
                packed[var] = [(str(i), v) for i, v in enumerate(callerlocal[var])]
            else:
                packed[var] = list(callerlocal[var].items())
                assert numkeys == len(packed[var]),\
                "A rule's input name must either all exist or not exist."

    with tempfile.NamedTemporaryFile(mode='wt') as tmpfile:
        json.dump(packed, tmpfile)
        tmpfile.flush()
        os.environ[PARAMETER_PASSING_ENVVAR] = tmpfile.name

        try:
            locals().update(callerglobal)
            locals().update(callerlocal)
            shell(_command)
        finally:
            del os.environ[PARAMETER_PASSING_ENVVAR]


class suffix_filter:

    def __init__(self, values):
        self.values = values

    def __getitem__(self, key):
        matches = [el for el in self.values if el.endswith(key)]
        if len(matches) > 1:
            raise ValueError("No single match found for {} in {}".format(key, self.values))
        elif len(matches) == 1:
            return matches[0]
        else:
            return ''


class temporary_file(str):

    def __new__(cls, suffix='', prefix='tmp', dir=None):
        tmpfile = tempfile.NamedTemporaryFile(suffix=suffix, prefix=prefix, dir=dir)

        self = str.__new__(cls, tmpfile.name)
        self.tmpfile = tmpfile
        return self

    def __del__(self):
        self.tmpfile.close()
        try:
            os.unlink(self.tmpfile.name)
        except:
            pass


def tmpfile(*args, **kwds):
    def finalize_temporary_file(_):
        return temporary_file(*args, **kwds)
    return finalize_temporary_file


try:
    import pushbullet

    PUSHBULLET_ENVVAR = 'SNAKEMAKE_PUSHBULLET_APIKEY'
    if PUSHBULLET_ENVVAR not in os.environ:
        raise ValueError("PushBullet API key is not defined in {}.".format(PUSHBULLET_ENVVAR))

    class PushBulletNotifier:

        def __init__(self, apikey, log_tail_length=30):
            self.apikey = apikey
            self.log_tail_length = log_tail_length
            self.pb = None

        def check_connection(self):
            if self.pb is None:
                self.pb = pushbullet.PushBullet(self.apikey)

        def message(self, title, msg, logfile=None):
            self.check_connection()

            self.pb.push_note(title, msg)
            if logfile is not None:
                logtail = open(logfile).readlines()[-self.log_tail_length:]
                self.pb.push_file(''.join(logtail), 'log.txt', file_type='text/plain')

    notify = PushBulletNotifier(os.environ[PUSHBULLET_ENVVAR])

except (ImportError, ValueError):
    class DummyNotifier:
        def __init__(self):
            pass
        def message(self, title, msg, logfile=None):
            pass
    notify = DummyNotifier()

# Call initializing functions
if is_snakemake_running():
    init_powersnake()
else:
    load_snakemake_params()

