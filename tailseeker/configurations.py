#!/usr/bin/env python3
#
# Copyright (c) 2015-6 Hyeshik Chang
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

    PATH_CONF_FILE = 'conf/paths.conf'

    def __init__(self, tailseeker_dir, settings_file):
        self.tailseeker_dir = tailseeker_dir
        self.confdata = self.load_config(settings_file)
        self.confdata = self.expand_sample_settings(self.confdata)
        self.paths = self.load_paths()

    def load_config(self, settings_file):
        usersettings = yaml.load(settings_file, Loader=yaml.FullLoader)
        if 'include' in usersettings:
            confdict = {}

            incfiles = usersettings['include']
            if isinstance(incfiles, str):
                incfiles = [incfiles]

            for predfile in incfiles:
                confpath = os.path.join(self.tailseeker_dir, 'conf', predfile)
                predconf = self.load_config(open(confpath))
                confdict = self.merge_configs(confdict, predconf)

            confdict = self.merge_configs(confdict, usersettings)
            return confdict
        else:
            return usersettings

    def merge_configs(self, conf1, conf2):
        merged = {}

        for conf1only in set(conf1) - set(conf2):
            merged[conf1only] = conf1[conf1only]
        for conf2only in set(conf2) - set(conf1):
            merged[conf2only] = conf2[conf2only]

        for shared in set(conf1) & set(conf2):
            if isinstance(conf1[shared], dict) != isinstance(conf1[shared], dict):
                raise ValueError('Layout is different between overriding configurations.')
            if isinstance(conf1[shared], dict):
                merged[shared] = self.merge_configs(conf1[shared], conf2[shared])
            else:
                merged[shared] = conf2[shared] # conf2 is overriding conf1.

        return merged

    def __getitem__(self, name):
        return self.confdata[name]

    def __contains__(self, name):
        return name in self.confdata

    def get(self, name, default=None):
        return self.confdata.get(name, default)

    def expand_sample_settings(self, node):
        predefined = {
            '_all': self.all_samples,
            '_exp': self.exp_samples,
            '_spk': self.spikein_samples,
        }

        finalized = {}
        for key, value in list(node.items()):
            if key not in predefined:
                finalized[key] = (node[key] if not isinstance(value, dict)
                                  else self.expand_sample_settings(value))
                continue

            for sample in predefined[key]:
                finalized[sample] = value

        node.clear()
        node.update(finalized)

        return finalized

    def load_paths(self):
        pathconf = os.path.join(self.tailseeker_dir, self.PATH_CONF_FILE)
        pathsettings = yaml.load(open(pathconf), Loader=yaml.FullLoader) if os.path.exists(pathconf) else {}
        if 'paths' in self.confdata:
            pathsettings.update(self.confdata['paths'])

        return pathsettings

    def export_paths(self, namespace, relative_to='/'):
        for progname, path in self.paths.items():
            resolved = os.path.abspath(os.path.join(relative_to, str(path)))
            namespace['{}_CMD'.format(progname.upper())] = resolved

    @property
    def all_samples(self):
        return sorted(map(str, list(self.exp_samples) +
                               list(self.spikein_samples)))

    @property
    def exp_samples(self):
        return sorted(self['experimental_samples'].keys())

    @property
    def spikein_samples(self):
        return sorted(self.get('spikein_samples', {}).keys())

    def get_sample_index(self, name):
        if name in self['experimental_samples']:
            return self['experimental_samples'][name]
        elif name in self.get('spikein_samples', {}):
            return self['spikein_samples'][name]
        else:
            raise KeyError("Sample {} not available in settings".format(repr(name)))


def scan_selectable_confs(tailseeker_dir):
    conffiles = glob.glob(os.path.join(tailseeker_dir, 'conf', '*.conf'))

    for conffilename in conffiles:
        confdata = yaml.load(open(conffilename), Loader=yaml.FullLoader)
        if 'name' in confdata:
            yield (confdata['name'], conffilename)

if __name__ == '__main__':
    print(list(scan_selectable_confs('..')))

