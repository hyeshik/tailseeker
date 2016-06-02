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


class CommandHandlers:

    COMMANDS = ['init', 'run', 'clean', 'clear']

    def __init__(self, tailseeker_dir):
        self.tailseeker_dir = tailseeker_dir
        self.tailseeker_main_file = os.path.join(tailseeker_dir, 'tailseeker', 'main.py')

    def proxy_to_snakemake(self, command, target):
        from snakemake import main

        sys.argv[0] = '{} {}'.format(sys.argv[0], command)
        sys.argv.insert(1, '-s')
        sys.argv.insert(2, self.tailseeker_main_file)
        if target:
            sys.argv.append(target)

        sys.exit(main())

    def check_configuration(self):
        if not os.path.exists('tailseeker.yaml'):
            print('ERROR: Configuration file "tailseeker.yaml" does '
                  'not exist in the current directory.', file=sys.stderr)
            sys.exit(1)

    def run(self):
        self.check_configuration()
        return self.proxy_to_snakemake('run', None)

    def clean(self):
        self.check_configuration()
        return self.proxy_to_snakemake('clean', 'clean')

    def clear(self):
        self.check_configuration()
        return self.proxy_to_snakemake('clear', 'clear')


def usage():
    print("""\
Tailseeker {version} - High-throughput measurement of poly(A) tails

Usage:  {command} <command> [options]

Commands:
  init      initiate and configure a project
  run       run an analysis workflow
  clean     clear intermediate files
  clear     clear all generated files
""".format(version=tailseeker.__version__, command=sys.argv[0]))


def main(tailseeker_dir):
    if len(sys.argv) < 2 or sys.argv[1] not in CommandHandlers.COMMANDS:
        usage()
        return

    os.environ['TAILSEEKER_DIR'] = tailseeker_dir

    command = sys.argv.pop(1)

    handlers = CommandHandlers(tailseeker_dir)
    getattr(handlers, command)()


if __name__ == '__main__':
    main()
