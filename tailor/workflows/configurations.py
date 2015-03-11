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

class Configurations:

    EXPANDABLE_CONF_KEYS = """
        spikein_training_length preamble_size preamble_sequence delimiter
        balance_check dupcheck_regions maximum_index_mismatches""".split()

    def __init__(self, settings_file):
        self.confdata = yaml.load(settings_file)
        self.expand_sample_settings()

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

