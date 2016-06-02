#!/usr/bin/env python3
#
# Copyright (c) 2016 Hyeshik Chang
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

from tailseeker.parsers import parse_sam, parse_taginfo, parse_fastq
from tailseeker.fileutils import MultiJoinIterator, TemporaryDirectory
from tailseeker.sequtils import reverse_complement_bytes, GiantFASTAFile
from tailseeker.parallel import open_tabix_parallel
from concurrent import futures
from scipy.stats import binom_test
import subprocess as sp
import sys
import os
import re

F_UNMAPPED          = 0x0004
F_REVERSE_STRAND    = 0x0010
F_SECOND_READ       = 0x0080

PAFLAG_CHIMERIC_TEMPLATE = 8192

SAMTOOLS_CMD = os.environ.get('SAMTOOLS_CMD', 'samtools')

cigar_pattern = re.compile(b'(\d+)([MIDNSHP])')
def calculate_cigar_length(cigar):
    lengths = [0, 0] # reference, read

    def found(match):
        size, action = match.groups()
        size = int(size)

        if action in b'MPDN': # was:MSPDN
            lengths[0] += size
        if action in b'MSHPI':
            lengths[1] += size

        return ''

    remaining = cigar_pattern.sub(found, cigar)
    if remaining != '':
        raise ValueError('Unknown pattern included: ' + remaining)

    return tuple(lengths)

def count_M_span(cigar):
    return sum(int(cnt) for cnt, cmd in cigar_pattern.findall(cigar) if cmd == 'M')

def cigar_reference_advance(cigar, pat=re.compile(b'\d+[MDN]')):
    return sum(int(tok[:-1]) for tok in pat.findall(cigar))

def test_terminal_match_seqs(refseq, readseq):
    if refseq == readseq:
        return 0

    for i, (refbase, rdbase) in enumerate(zip(refseq, readseq)):
        if refbase != rdbase:
            return len(refseq) - i
    return 0

MAPPED_READ_MASK = F_UNMAPPED | F_SECOND_READ
is_r5_mapped = lambda x: (x.flag & MAPPED_READ_MASK) == 0
is_mapped = lambda x: (x.flag & F_UNMAPPED) != F_UNMAPPED
is_r3 = lambda x: (x.flag & F_SECOND_READ) == F_SECOND_READ

def reevaluate_terminal_additions(samrows, taginforows, refgenome,
                                  maximum_fragment_length, terminal_recheck_width):
    rflags = 0
    samrows = list(samrows)
    taginforows = list(taginforows)
    if not taginforows:
        return

    taginfo = taginforows[0]

    read3rows = list(filter(is_r3, samrows))
    if not read3rows:
        raise ValueError("Unmatched read pair {} found.".format(samrows[0].qname.decode()))

    read3rows_mapped = list(filter(is_mapped, read3rows))

    # remove potential artifacts from RNA-RNA ligations
    if read3rows_mapped and maximum_fragment_length is not None:
        read5rows = list(filter(is_r5_mapped, samrows))

        valid_read3rows = [
            r3row for r3row in read3rows_mapped
            if any((r5row.rname == r3row.rname and
                        abs(r5row.pos - r3row.pos) < maximum_fragment_length)
                    for r5row in read5rows)
        ]

        if read3rows_mapped and not valid_read3rows:
            # likely a chimeric clone
            rflags |= PAFLAG_CHIMERIC_TEMPLATE
        #else:
            # This filter appears to produce more false positives than expected
            # purifications.
            # read3rows_mapped = valid_read3rows

    if not read3rows_mapped: # all read 3 alignments were unmapped or filtered out
        fullmodlen = taginfo.polyA + len(taginfo.mods)
        r3seq = read3rows[0].seq
        if (read3rows[0].flag & F_REVERSE_STRAND) == 0:
            r3seq = reverse_complement_bytes(r3seq)
        modseq = r3seq[-fullmodlen:].decode() if fullmodlen > 0 else ''

        return taginfo, '-', modseq, rflags

    read3rows_mapped.sort(key=lambda x: (-count_M_span(x.cigar), -x.mapq))
    bestmap = read3rows_mapped[0]
    additional_clipping = 0

    cigartokens = cigar_pattern.findall(bestmap.cigar)
    if bestmap.flag & F_REVERSE_STRAND:
        # read 2 is reverse mapped -> RNA is same to + strand
        modlen = int(cigartokens[-1][0]) if cigartokens[-1][1] == b'S' else 0

        # realign the terminal M ops.
        mapped_end = bestmap.pos - 1 + cigar_reference_advance(bestmap.cigar)
        matchseq = bestmap.seq[:-modlen] if modlen else bestmap.seq
        refseq = refgenome.get(bestmap.rname.decode(),
                               mapped_end - terminal_recheck_width,
                               mapped_end).upper()
        if (len(refseq) == terminal_recheck_width and
                len(matchseq) >= terminal_recheck_width):
            additional_clipping = test_terminal_match_seqs(refseq.encode(),
                                                matchseq[-terminal_recheck_width:])
            if additional_clipping > 0:
                modlen += additional_clipping

        # common stuff for + strand
        modseq = bestmap.seq[-modlen:] if modlen else b''
        strand = b'+'
    else:
        # read 2 is forward mapped -> RNA is same to - strand
        modlen = int(cigartokens[0][0]) if cigartokens[0][1] == b'S' else 0

        # realign the terminal M ops.
        mapped_end = bestmap.pos - 1
        matchseq = reverse_complement_bytes(bestmap.seq[modlen:] if modlen else bestmap.seq)
        refseq = refgenome.get(bestmap.rname.decode(), mapped_end,
                               mapped_end + terminal_recheck_width, '-').upper()
        if (len(refseq) == terminal_recheck_width and
                len(matchseq) >= terminal_recheck_width):
            additional_clipping = test_terminal_match_seqs(refseq.encode(),
                                                matchseq[-terminal_recheck_width:])
            if additional_clipping > 0:
                modlen += additional_clipping

        # common stuff for - strand
        modseq = reverse_complement_bytes(bestmap.seq[:modlen]) if modlen else b''
        strand = b'-'

    cigarpart = strand + bestmap.cigar

    return taginfo, cigarpart.decode(), modseq.decode(), rflags


def calculate_required_A_in_polyA(prob, pcutoff, maximum_length):
    required_A = []

    for i in range(maximum_length + 1):
        best = 0
        for j in range(i):
            if binom_test(j, i, prob, alternative='less') < pcutoff:
                best = j
            else:
                break
        required_A.append(i - best)

    return required_A


def process_tile(options, tile, taginfo_input, output):
    readerproc = sp.Popen([SAMTOOLS_CMD, 'view', '-r', tile, options.aln], stdout=sp.PIPE)

    # open the reference database
    refgenome = GiantFASTAFile(options.refseq_fasta)

    # parsing iterators
    samit = parse_sam(readerproc.stdout)
    taginfoit = parse_taginfo(taginfo_input())

    # key functions for joiner
    def parse_readid_from_sam(samrow):
        qtokens = samrow.qname.split(b':', 2)
        return (qtokens[0], int(qtokens[1]))
    taginfokey = lambda x: (x.tile, x.cluster)

    joined_it = MultiJoinIterator([samit, taginfoit], [parse_readid_from_sam, taginfokey])
    re_termmod = re.compile('(T+|G+|C+)?$')
    re_termA = re.compile('A*$')

    polya_noreev_threshold = options.maxpolyareev
    required_A = calculate_required_A_in_polyA(options.reev_cprob, 0.95, 100)
    rescue_threshold = options.rescue_threshold

    outfile = open(output, 'w')

    for (tile, cluster), samrows, taginforows in joined_it:
        samrows = list(samrows)
        if len(samrows) == 0:
            continue # may filtered out by the contaminant filter.

        taginfo, cigar, nontmplmods, rflags = reevaluate_terminal_additions(
                samrows, taginforows, refgenome, options.fragsize,
                options.termaln_check)

        if taginfo.polyA > polya_noreev_threshold:
            polyA_len = taginfo.polyA
        else:
            tailbody = re_termmod.sub('', nontmplmods)
            if not tailbody:
                polyA_len = 0
            elif tailbody.count('A') >= required_A[len(tailbody)]:
                polyA_len = tailbody.count('A')
            else:
                termA_len = len(re_termA.findall(tailbody)[0])
                if termA_len >= rescue_threshold:
                    polyA_len = termA_len
                else:
                    polyA_len = 0

        reev_mods = {'A': polyA_len, 'G': 0, 'C': 0, 'T': 0}

        if nontmplmods:
            if polyA_len >= 3:
                termmods = re_termmod.findall(nontmplmods)
                if termmods[0]:
                    reev_mods[termmods[0][0]] = len(termmods[0])
            elif len(set(nontmplmods)) == 1 and nontmplmods[0] != 'A':
                reev_mods[nontmplmods[0]] = len(nontmplmods)

        print(tile.decode(), cluster, taginfo.pflags | rflags, taginfo.clones,
              reev_mods['A'] if taginfo.polyA >= 0 else -1,
              reev_mods['T'], reev_mods['G'], reev_mods['C'],
              nontmplmods, sep='\t', file=outfile)

    return output


def process(options):
    taginfo_inputs = open_tabix_parallel(options.taginfo, named=True)

    with futures.ProcessPoolExecutor(options.parallel) as executor, \
            TemporaryDirectory(asobj=True) as workdir:
        jobs = []

        for tile, taginfo_open in sorted(taginfo_inputs.items()):
            joboutput = os.path.join(workdir.path, tile)
            job = executor.submit(process_tile, options, tile, taginfo_open, joboutput)
            jobs.append(job)

        for j in jobs:
            joboutput = j.result()
            sp.check_call(['cat', joboutput], stdout=sys.stdout)
            os.unlink(joboutput)


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Generate a table for '
                                                 'non-templated 3\' additions')
    parser.add_argument('--parallel', dest='parallel', metavar='NUM', type=int,
                        default=1, help='Number of parallel processes')
    parser.add_argument('--taginfo', dest='taginfo', metavar='FILE', type=str,
                        required=True, help='Path to a taginfo file')
    parser.add_argument('--alignment', dest='aln', metavar='FILE', type=str,
                        required=True, help='Path to a paired alignment file (bam)')
    parser.add_argument('--reference-seq', dest='refseq_fasta', metavar='FILE', type=str,
                        required=True, help='Path to a reference genome FASTA file')
    parser.add_argument('--max-fragment-size', dest='fragsize', metavar='SIZE',
                        type=int, default=None,
                        help='Maximum distance between two reads to suppress RNA-RNA '
                             'in vitro fusion artifacts.')
    parser.add_argument('--terminal-mismatch-recheck', dest='termaln_check', metavar='WIDTH',
                        type=int, default=2,
                        help="Width of regions in M ops that will be rechecked if "
                             "containing mismatches and converted to soft clips.")
    parser.add_argument('--max-polya-reevaluation', dest='maxpolyareev', metavar='LENGTH',
                        type=int, default=7,
                        help='Maximum length of poly(A) to reevaluate '
                             'with the reference sequence')
    parser.add_argument('--polya-reev-contamination-prob', dest='reev_cprob', metavar='PROB',
                        type=float, default=.2,
                        help='Probability of contaminant in poly(A) tails to be re-evaluated')
    parser.add_argument('--contaminated-polya-rescue-threshold', dest='rescue_threshold',
                        metavar='LENGTH', type=int, default=3,
                        help='Length of A-stretches at the 3\'-end that enforces A length '
                             'measurement even though the untemplated tailing as whole is '
                             'contaminanted by other type of tails.')
    options = parser.parse_args()

    return options

if __name__ == '__main__':
    options = parse_arguments()
    process(options)
