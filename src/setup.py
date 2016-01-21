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

import os
import subprocess as sp
from distutils.core import setup, Extension

if 'TAILSEEKER_DIR' in os.environ:
    TAILSEEKER_DIR = os.environ['TAILSEEKER_DIR']
else:
    TAILSEEKER_DIR = os.path.abspath('..')

PKGCONFIG_CMD = os.environ.get('PKGCONFIG_CMD', 'pkg-config')

# This function implementation was borrowed from George Snyder.
# - http://code.activestate.com/recipes/502261-python-distutils-pkg-config/#c2
def pkgconfig(package, **kw):
    flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
    for token in sp.check_output("{} --libs --cflags {}".format(PKGCONFIG_CMD, package),
                                 shell=True).decode('ascii').split():
        if token[:2] in flag_map:
            kw.setdefault(flag_map.get(token[:2]), []).append(token[2:])
        else:
            kw.setdefault('extra_compile_args', []).append(token)
    return kw

def fix_tailseeker_path():
    content = open('tailseeker.in').read()
    content_fixed = content.replace('%%TAILSEEKER_DIR%%', TAILSEEKER_DIR)
    open('tailseeker', 'w').write(content_fixed)

fix_tailseeker_path()

srcfile = 'pythonext/tailseqext.c'
modname = 'tailseqext'
htslib_flags = pkgconfig('htslib')

setup(name = modname,
      version = '0.2',
      description = 'Internal extension module for tailseeker',
      author = 'Hyeshik Chang',
      author_email = 'hyeshik@snu.ac.kr',
      url = 'http://highthroughput.org',
      license = 'MIT',
      ext_modules = [Extension(modname, [srcfile], **htslib_flags)],
      scripts = ['tailseeker']
)
