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

import sys
import os
import shutil


def setdefault(name, value):
    if name not in globals():
        globals()[name] = value

def make_scratch_dir(prefix):
    dirprefix = os.path.join(SCRATCHDIR, prefix)
    if not os.path.isdir(dirprefix):
        os.makedirs(dirprefix)
    return dirprefix

def create_scratch_link():
    linkpath = os.path.abspath(os.path.join(WRKDIR, 'scratch'))
    origpath = os.path.abspath(SCRATCHDIR)
    if not os.path.exists(linkpath) and linkpath != origpath:
        if not os.path.isdir(origpath):
            os.makedirs(origpath)
        os.symlink(origpath, linkpath)

def init_default_dirs():
    setdefault('WRKDIR', os.getcwd())

    if 'TOPDIR' not in globals():
        if 'TAILSEQDIR' in os.environ:
            TOPDIR = os.environ['TAILSEQDIR']
        else:
            TOPDIR = find_topdir(WRKDIR)
    RESOURCEDIR = os.path.join(TOPDIR, 'resources')
    LOCALRESOURCEDIR = os.path.join(TOPDIR, 'local-resources')
    TOOLSDIR = os.path.join(TOPDIR, 'tools')
    setdefault('SCRATCHDIR', os.path.join(WRKDIR, 'scratch'))
    create_scratch_link()

class SuffixFilter:
    def __init__(self, values):
        self.values = values

    def __getitem__(self, key):
        matches = [el for el in self.values if el.endswith(key)]
        if len(matches) != 1:
            raise ValueError("No single match found for {} in {}".format(key, self.values))
        return matches[0]



#========== Initializations ============
sys.path.append(TAILOR_DIR)

# Load and parse configuration settings.
SETTINGS_FILE = os.path.abspath('tailorconf.yaml')

from tailor import configurations
CONF = configurations.Configurations(TAILOR_DIR, open(SETTINGS_FILE))

# Verify directories and links.
WRKDIR = os.getcwd()
RESOURCESDIR = os.path.join(TAILOR_DIR, 'resources')
BINDIR = os.path.join(TAILOR_DIR, 'bin')
SCRATCHDIR = (CONF['scratch_dir'] if 'scratch_dir' in CONF
                                  else os.path.join(TAILOR_DIR, 'scratch'))
SCRIPTSDIR = os.path.join(TAILOR_DIR, 'scripts')
CONF.export_paths(globals())
create_scratch_link()

if not CONF['clean_intermediate_files']:
    nonfinal = temp = lambda arg: arg
elif CONF['clean_intermediate_files'] == 'some':
    nonfinal = lambda arg: arg
elif CONF['clean_intermediate_files'] == 'immediately':
    nonfinal = temp
else:
    raise ValueError('Setting clean_intermediate_files should be one among immediately, some, '
                     'or no.')


# XXX: move these parameters somewhere out
SIGNAL_STABILIZER_TARGET_RANGE = 1, 50 # in 1-based inclusive coordinate
SIGNAL_STABILIZER_REFERENCE_RANGE = 51, 55 # 1-based
SIGNAL_STABILIZER_POLYA_DETECTION_RANGE = 51, 100 # 1-based

PASIGNAL_CLIP_MIN = -0.3
PASIGNAL_CLIP_MAX = 2.5

MAXIMUM_DELIMITER_MISALIGNMENT = 1
# XXX

# Predefined constants
inf = float('inf')
nan = float('nan')

# Commands needs to be run with bash with these options to terminate on errors correctly.
shell.executable(os.popen('which bash').read().strip()) # pipefail is supported by bash only.
shell.prefix(('set -e; set -o pipefail; '
              'export PYTHONPATH="{PYTHONPATH}" '
                     'BGZIP_CMD="{HTSLIB_BINDIR}/bgzip" '
                     'TABIX_CMD="{HTSLIB_BINDIR}/tabix" '
                     'TAILSEQ_SCRATCH_DIR="{SCRATCHDIR}" '
                     + CONF.get('envvars', '') + '; ').format(
                PYTHONPATH=TAILOR_DIR, SCRATCHDIR=SCRATCHDIR, HTSLIB_BINDIR=HTSLIB_BINDIR))

