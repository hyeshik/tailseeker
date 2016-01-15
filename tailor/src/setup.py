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
from distutils.core import setup, Extension

if 'TAILOR_DIR' in os.environ:
    TAILOR_DIR = os.environ['TAILOR_DIR']
else:
    TAILOR_DIR = os.path.abspath('..')


def fix_tailor_path():
    content = open('tailseeker.in').read()
    content_fixed = content.replace('%%TAILOR_DIR%%', TAILOR_DIR)
    open('tailseeker', 'w').write(content_fixed)

fix_tailor_path()

srcfile = 'pythonext/tailseqext.c'
modname = 'tailseqext'

setup(name = modname,
      version = '0.2',
      description = 'extension modules for TAIL-seq',
      author = 'Hyeshik Chang',
      author_email = 'hyeshik@snu.ac.kr',
      url = 'http://highthroughput.org',
      license = 'MIT',
      ext_modules = [Extension(modname, [srcfile])],
      scripts = ['tailseeker']
)
