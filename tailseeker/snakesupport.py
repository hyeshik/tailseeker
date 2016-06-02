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
from collections import defaultdict
from snakemake.io import limit


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

all_intermediate_files = set()
def lazy_clearing(arg):
    all_intermediate_files.add(arg)
    return arg

def immediate_clearing(arg, sntemp=temp):
    all_intermediate_files.add(arg)
    return sntemp(arg)


#========== Initializations ============
from tailseeker.powersnake import *

# Load and parse configuration settings.
SETTINGS_FILE = os.path.abspath('tailseeker.yaml')

from tailseeker import configurations
CONF = configurations.Configurations(TAILSEEKER_DIR, open(SETTINGS_FILE))

# Verify directories and links.
WRKDIR = os.getcwd()
RESOURCESDIR = os.path.join(TAILSEEKER_DIR, 'resources')
BINDIR = os.path.join(TAILSEEKER_DIR, 'bin')
SCRATCHDIR = (CONF['scratch_dir'] if 'scratch_dir' in CONF
                                  else os.path.join(TAILSEEKER_DIR, 'scratch'))
SCRIPTSDIR = os.path.join(TAILSEEKER_DIR, 'scripts')
CONF.export_paths(globals())
create_scratch_link()

if not CONF['clean_intermediate_files']:
    nonfinal = temp = lazy_clearing
elif CONF['clean_intermediate_files'] == 'some':
    nonfinal, temp = lazy_clearing, immediate_clearing
elif CONF['clean_intermediate_files'] == 'immediately':
    nonfinal = temp = immediate_clearing
else:
    raise ValueError('Setting clean_intermediate_files should be one among immediately, some, '
                     'or no.')


# Predefined constants
inf = float('inf')
nan = float('nan')

# Commands needs to be run with bash with these options to terminate on errors correctly.
shell.executable(BASH_CMD) # pipefail is supported by bash only.
shell.prefix(('set -e; set -o pipefail; '
              'export PYTHONPATH="{PYTHONPATH}" LC_ALL=C '
                     'BGZIP_CMD="{BGZIP_CMD}" '
                     'TABIX_CMD="{TABIX_CMD}" '
                     'TAILSEQ_SCRATCH_DIR="{SCRATCHDIR}" '
                     + CONF.get('envvars', '') + '; ').format(
                PYTHONPATH=TAILSEEKER_DIR, SCRATCHDIR=SCRATCHDIR, BGZIP_CMD=BGZIP_CMD,
                TABIX_CMD=TABIX_CMD))

