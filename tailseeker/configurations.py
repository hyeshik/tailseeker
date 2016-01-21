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

import yaml
import os
import glob


class Configurations:

    EXPANDABLE_CONF_KEYS = """
        spikein_training_length umi_length umi_fixed_sequence delimiter
        balance_check dupcheck_regions maximum_index_mismatches
        species""".split()

    PATH_CONF_FILE = 'conf/paths.conf'

    def __init__(self, tailor_dir, settings_file):
        self.tailor_dir = tailor_dir
        self.confdata = self.load_config(settings_file)
        self.expand_sample_settings()

    def load_config(self, settings_file):
        usersettings = yaml.load(settings_file)
        if 'include' in usersettings:
            confdict = {}

            incfiles = usersettings['include']
            if isinstance(incfiles, str):
                incfiles = [incfiles]

            for predfile in incfiles:
                confpath = os.path.join(self.tailor_dir, 'conf', predfile)
                predconf = self.load_config(open(confpath))
                confdict.update(predconf)

            confdict.update(usersettings)
            return confdict
        else:
            return usersettings

    def __getitem__(self, name):
        return self.confdata[name]

    def __contains__(self, name):
        return name in self.confdata

    def get(self, name, default=None):
        return self.confdata.get(name, default)

    def expand_sample_settings(self):
        predefined = {
            '_all': self.all_samples,
            '_exp': self.exp_samples,
            '_spk': self.spikein_samples,
        }

        for key in self.EXPANDABLE_CONF_KEYS:
            if key not in self.confdata:
                continue

            user_specified = self.confdata[key]
            finalized = {}

            for k, d in user_specified.items():
                if k in predefined:
                    for s in predefined[k]:
                        finalized[s] = d
                else:
                    finalized[k] = d

            user_specified.clear()
            user_specified.update(finalized)

    def export_paths(self, namespace):
        pathconf = os.path.join(self.tailor_dir, self.PATH_CONF_FILE)
        pathsettings = yaml.load(open(pathconf)) if os.path.exists(pathconf) else {}
        if 'paths' in self.confdata:
            pathsettings.update(self.confdata['paths'])

        for progname, path in pathsettings.items():
            namespace['{}_CMD'.format(progname.upper())] = path

    @property
    def all_samples(self):
        return sorted(map(str, list(self['experimental_samples'].keys()) +
                               list(self['spikein_samples'].keys())))

    @property
    def exp_samples(self):
        return sorted(self['experimental_samples'].keys())

    @property
    def spikein_samples(self):
        return sorted(self['spikein_samples'].keys())

    def get_sample_index(self, name):
        if name in self['experimental_samples']:
            return self['experimental_samples'][name]
        elif name in self['spikein_samples']:
            return self['spikein_samples'][name]
        else:
            raise KeyError("Sample {} not available in settings".format(repr(name)))


def scan_selectable_confs(tailor_dir):
    conffiles = glob.glob(os.path.join(tailor_dir, 'conf', '*.conf'))

    for conffilename in conffiles:
        confdata = yaml.load(open(conffilename))
        if 'name' in confdata:
            yield (confdata['name'], conffilename)

if __name__ == '__main__':
    print(list(scan_selectable_confs('..')))

