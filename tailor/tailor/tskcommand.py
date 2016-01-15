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
import tailor


class CommandHandlers:

    COMMANDS = ['init', 'run', 'clean']

    def __init__(self, tailor_dir):
        self.tailor_dir = tailor_dir
        self.tailor_main_file = os.path.join(tailor_dir, 'tailor', 'main.py')

    def run(self):
        from snakemake import main

        sys.argv[0] = sys.argv[0] + ' run'
        sys.argv.insert(1, '-s')
        sys.argv.insert(2, self.tailor_main_file)
        sys.exit(main())

    def clean(self):
        from snakemake import main

        sys.argv[0] = sys.argv[0] + ' clean'
        sys.argv.insert(1, '-s')
        sys.argv.insert(2, self.tailor_main_file)
        sys.argv.append('clean')

        sys.exit(main())


def usage():
    print("""\
Tailseeker {version} - High-throughput measurement of poly(A) tails

Usage:  {command} <command> [options]

Commands:
  init      initiate and configure a project
  run       run an analysis workflow
  clean     clear intermediate files
""".format(version=tailor.__version__, command=sys.argv[0]))


def main(tailor_dir):
    if len(sys.argv) < 2 or sys.argv[1] not in CommandHandlers.COMMANDS:
        usage()
        return

    os.environ['TAILOR_DIR'] = tailor_dir

    command = sys.argv.pop(1)

    handlers = CommandHandlers(tailor_dir)
    getattr(handlers, command)()


if __name__ == '__main__':
    main()
