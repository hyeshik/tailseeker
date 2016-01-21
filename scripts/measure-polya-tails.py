#!/usr/bin/env python3
#
# Copyright (c) 2013 Institute for Basic Science
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
# - Hyeshik Chang <hyeshik@snu.ac.kr>
#

from tailor.parsers import parse_sqi, parse_pascore
from tailor.parallel import open_tabix_parallel
from tailor.fileutils import (
    ParallelMatchingReader, TemporaryDirectory, MultiJoinIterator)
from tailseqext import SimpleBGZFWriter
from concurrent import futures
import ghmm
import csv
import re
import numpy as np
from collections import Counter
from itertools import takewhile
import pickle


class HMMPolyARuler(object):

    SENTINEL = 1000000
    POLYA_STATES = (0, 1)

    def __init__(self, modelfile):
        modeldata = pickle.load(open(modelfile, 'rb'))
        self.domain = ghmm.Float()
        self.model = ghmm.HMMFromMatrices(self.domain,
                                          ghmm.GaussianMixtureDistribution(self.domain),
                                          *modeldata['parameters'])
        self.vrange = modeldata['vrange']

    def rule(self, pascore):
        seq = ghmm.EmissionSequence(self.domain,
                                    [self.SENTINEL] +
                                    list(map(float, pascore.clip(*self.vrange))))
        vitstates = self.model.viterbi(seq)[0]
        return self.determine_length(vitstates[1:])

    def determine_length(self, states, polya_states=POLYA_STATES):
        return len(list(takewhile(lambda x: x in polya_states, states)))


def rule_polya_naive(r2seq, pat=re.compile('^T*')):
    return len(pat.findall(r2seq)[0])


def print_debug_spot(sqispot, paspot, seq_start_pa, seq_end_pa, hmm_length):
    seq = sqispot.seq[sqispot.istart:]
    qual = sqispot.qual[sqispot.istart:]
    pa = paspot.pascore[-len(seq):]

    import csv
    out = csv.writer(open('debug.csv', 'w'))
    out.writerow(['seq'] + list(seq[seq_start_pa:]))
    out.writerow(['qual'] + list(qual[seq_start_pa:]))
    out.writerow(['pa'] + list(pa[seq_start_pa:]))

    read2umipart = slice(57, 72)
    normbasis = sqispot.intensity[read2umipart].mean(axis=0).clip(1)
    normintensity = sqispot.intensity[sqispot.istart + seq_start_pa:].clip(1) / normbasis

    for chan, intensity in zip('ACGT', sqispot.intensity[sqispot.istart + seq_start_pa:].transpose()):
        out.writerow([chan] + list(intensity))

    for chan, intensity in zip('ACGT', normintensity.transpose()):
        out.writerow(['norm ' + chan] + list(intensity))

    for chan, intensity in zip('ACGT', sqispot.intensity[read2umipart].transpose()):
        out.writerow(['umi ' + chan] + list(intensity))

    # rescaling method
    print(sqispot.intensity[read2umipart])
    scorebins = dict((nt, [[], []]) for nt in 'ACGTN')
    for base, intensity in zip(sqispot.seq[read2umipart], sqispot.intensity[read2umipart]):
        for ib, ivalue in zip('ACGT', intensity):
            scorebins[ib][base == ib].append(ivalue)

    print(scorebins)
    origintensity = sqispot.intensity[sqispot.istart + seq_start_pa:].transpose()
    for i, nt in enumerate('ACGT'):
        neg = np.log2(np.array(scorebins[nt][0]).clip(1)).mean()
        pos = np.log2(np.array(scorebins[nt][1]).clip(1)).mean()

        print(nt, neg, pos)
        rescaled = (np.log2(origintensity[i].clip(1)) - neg) / (pos - neg)
        out.writerow(['rescaled ' + nt] + list(rescaled))


def run_subjob(sqiopener, paopener, options, joboutput):
  try:
    sqi_in = parse_sqi(sqiopener())
    pa_in = parse_pascore(paopener())
    common_key = lambda s: (s.tile, s.cluster)

    preader = MultiJoinIterator([sqi_in, pa_in], [common_key, common_key])
    writer = SimpleBGZFWriter(joboutput)
    hmmruler = HMMPolyARuler(options.model)

    for (tile, cluster), sqispot, paspot in preader:
        sqispot, paspot = next(sqispot), next(paspot, None)

        seq = sqispot.seq[sqispot.istart:]

        seq_start_pa = paspot.endmod_len if paspot is not None else 0
        seqbased_length = paspot.seqbased_polya_len if paspot is not None else 0
        hmm_length = hmmruler.rule(paspot.pascore) if paspot is not None else 0
        naive_length = rule_polya_naive(seq)

        if paspot is not None:
            # HMM often mis-calls short poly(A) tails as zero-length.
            final_length = hmm_length if hmm_length > 4 else seqbased_length
        else:
            final_length = seqbased_length

        outline = '\t'.join(map(str, (tile, cluster, seq_start_pa, final_length,
                                      seqbased_length, hmm_length, naive_length))) + '\n'
        writer(outline)

    writer.close()

  except:
    import traceback
    traceback.print_exc()
    raise


def terminated_iterator():
    return
    yield None


def run(options):
    sqi_openers = open_tabix_parallel(options.input_sqi, named=True)
    pa_openers = open_tabix_parallel(options.input_pa, named=True)

    tiles = sorted(sqi_openers)
    sqi_openers = [sqi_openers[k] for k in tiles]
    pa_openers = [(pa_openers[k] if k in pa_openers else terminated_iterator) for k in tiles]

    model = pickle.load(open(options.model, 'rb'))

    with futures.ProcessPoolExecutor(options.parallel) as executor, \
            TemporaryDirectory(asobj=True) as workdir:
        jobs = []

        for inputno, (sqiopen, paopen) in enumerate(zip(sqi_openers, pa_openers)):
            joboutput = '{}/{:08x}'.format(workdir.path, inputno)
            job = executor.submit(run_subjob, sqiopen, paopen, options, joboutput)
            jobs.append(job)

        for j in jobs:
            j.result()

        workdir.merge_bgzf(options.output)



def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Locate Poly(A) tails in sequence reads and '
                                                'PA scores')
    parser.add_argument('--input-sqi', dest='input_sqi', metavar='FILE', required=True,
                        type=str, help='Path to a file containing sequence and intensity')
    parser.add_argument('--input-pa', dest='input_pa', metavar='FILE', required=True,
                        type=str, help='Path to a file containing PA scores')
    parser.add_argument('--model', dest='model', metavar='FILE', required=True,
                        type=str, help='Path to a file containing poly(A) tail GMHMM model')
    parser.add_argument('--score-T', dest='score_t', metavar='N',
                        default=1, type=int, help='Alignment score for T (for v1 only)')
    parser.add_argument('--score-N', dest='score_n', metavar='N',
                        default=-5, type=int, help='Alignment score for N (for v1 only)')
    parser.add_argument('--score-ACG', dest='score_acg', metavar='N',
                        default=-10, type=int, help='Alignment score for A, C or G '
                                                    '(for v1 only)')
    parser.add_argument('--max-seq-counting', dest='max_seq_counting', metavar='N',
                        default=8, type=int,
                        help='Maximum poly(A) length to count based on sequences '
                             '(without HMM) (for v1 only)')
    parser.add_argument('--max-terminal-modification', dest='max_term_mod', metavar='N',
                        default=30, type=int,
                        help='Maximum length of 3\' terminal modification after '
                             'poly(A) tails (for v1 only)')
    parser.add_argument('--parallel', dest='parallel', metavar='N',
                        default=8, type=int, help='Number of parallel processes')
    parser.add_argument('--output', dest='output', metavar='FILE', type=str,
                        required=True, help='Path to output file (bgz format)')

    options = parser.parse_args()

    return options


if __name__ == '__main__':
    options = parse_arguments()
    run(options)

