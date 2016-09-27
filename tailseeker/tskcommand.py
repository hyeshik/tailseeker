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
import sys
import tailseeker


def check_configuration():
    if not os.path.exists('tailseeker.yaml'):
        print('ERROR: Configuration file "tailseeker.yaml" does '
              'not exist in the current directory.', file=sys.stderr)
        sys.exit(1)

def proxy_to_snakemake(tailseeker_dir):
    from snakemake import main

    tailseeker_main_file = os.path.join(tailseeker_dir, 'tailseeker', 'main.py')

    if '-s' not in sys.argv and '--snakefile' not in sys.argv:
        sys.argv.insert(1, '-s')
        sys.argv.insert(2, tailseeker_main_file)

    sys.exit(main())

def show_banner():
    print("""\
Tailseeker {version} - High-throughput measurement of poly(A) tails
""".format(version=tailseeker.__version__))

def main(tailseeker_dir):
    show_banner()

    os.environ['TAILSEEKER_DIR'] = tailseeker_dir

    check_configuration()
    proxy_to_snakemake(tailseeker_dir)


if __name__ == '__main__':
    main()
