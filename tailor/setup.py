#!/usr/bin/env python
import sys
from distutils.core import setup, Extension

TAILSEQ_MEGADEPS = '/home/tsuser/ts-deps'
NUMPY_INCLUDE = '{prefix}/lib/{ver}/site-packages/numpy/core/include'.format(
                       prefix=TAILSEQ_MEGADEPS, ver='python{}.{}'.format(*sys.version_info[:2]))

setup (name = 'tailseqext',
       version = '0.1',
       description = 'extension modules for tail-seq',
       author = 'Hyeshik Chang',
       author_email = 'hyeshik@snu.ac.kr',
       url = 'http://www.narrykim.org',
       license = 'MIT',
       ext_modules = [Extension('tailseqext', ['tailseqext.c'])],
       include_dirs = [NUMPY_INCLUDE]
)
