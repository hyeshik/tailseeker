#!/usr/bin/env python3
#
# Copyright (c) 2013-2015 Institute for Basic Science
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

from tailseeker.fileutils import TemporaryDirectory, MultiJoinIterator
from tailseeker.sequencers import decode_gcid, tile_range_cond
import os
import sys
import csv
import tables
import subprocess as sp
import pandas as pd
from itertools import groupby
from concurrent import futures
from functools import partial
from collections import defaultdict
from random import randrange, seed
from io import TextIOWrapper


class DuplicateCount(tables.IsDescription):
    gcid = tables.Int64Col(pos=0)
    count = tables.Int64Col(pos=1)


class SequenceFilter:

    def __init__(self, name):
        self.passed = self.failed = 0
        self.name = name

    def get_stats(self):
        return self.name, self.passed, self.failed

    def __call__(self, seq):
        passed = self.perform(seq)

        if passed:
            self.passed += 1
        else:
            self.failed += 1

        return passed


class BalanceFilter(SequenceFilter):

    def __init__(self, options):
        if options.balanceregion is None:
            raise ValueError("Options for BalanceFilter is not available.")

        super().__init__('umi_balance')
        balance_start1, balance_end = map(int, options.balanceregion.split(':'))
        self.balanceslice = slice(balance_start1 - 1, balance_end)
        self.balancemin = options.balancemin

    def perform(self, seq):
        bseq = seq[self.balanceslice]
        return min(bseq.count(c) for c in b'ACGT') >= self.balancemin


class FixedUMIFilter(SequenceFilter):

    def __init__(self, options):
        if options.umipos is None:
            raise ValueError("Options for BalanceFilter is not available.")

        super().__init__('fixed_umi')
        self.umipos, self.umiseq = options.umipos - 1, options.umiseq.encode()
        self.umislice = slice(self.umipos, self.umipos + len(self.umiseq))
        self.umimismatch = options.umimismatch
        self.passed = self.failed = 0

    @staticmethod
    def string_dist(s1, s2):
        return sum(a != b for a, b in zip(s1, s2))

    def perform(self, seq):
        return self.string_dist(seq[self.umislice], self.umiseq) <= self.umimismatch


def generate_id_sequences(input, samplename, num, slices, tmpdir):
    outputfile = os.path.join(tmpdir, str(num))
    output = open(outputfile, 'wb')

    with tables.open_file(input, 'r') as h5file:
        seqquals = h5file.get_node('/primary_sqi/{}/seqqual'.format(samplename))
        for spot in seqquals:
            seq = spot['seq']
            idseq = b''.join(seq[s] for s in slices)
            output.write(idseq)
            output.write('\t{}\n'.format(spot['gcid']).encode('ascii'))

    output.close()

    return outputfile


def sort_and_uniquify_duplicates(output, samplename, idfiles, parallel, debug_output=None,
                                 clonelist_output=None):
    h5 = tables.open_file(output, 'w')
    group = h5.create_group('/', 'uniqids', 'Deduplicated ID list of unique clusters')
    uniqids_out = h5.create_table(group, samplename + '_working', DuplicateCount,
                                  'Deduplicated ID list for sample {}'.format(samplename))

    soproc = sp.Popen(['sort', '--parallel={}'.format(parallel), '-k1,1'] + idfiles,
                      stdout=sp.PIPE)
    sofeedout = (line.split() for line in soproc.stdout)
    grpstats = defaultdict(int)

    seed(8809)

    for seq, grp in groupby(sofeedout, key=lambda x: x[0]):
        grp = [int(gcid.decode()) for _, gcid in grp]

        duplicity = len(grp)
        grpstats[duplicity] += 1

        if duplicity >= 2:
            survivor = randrange(duplicity)
            surviving_cluster = grp.pop(survivor)

            if debug_output is not None and len(grp) >= 10:
                print('{}\t{}\t{}'.format(duplicity, seq.decode(), surviving_cluster),
                      file=debug_output)
        else:
            surviving_cluster = grp.pop()

        if clonelist_output is not None:
            print(surviving_cluster, ' '.join(map(str, grp)), sep='\t',
                  file=clonelist_output)

        uniqids_out.append([(surviving_cluster, duplicity)])

    uniqids_out.cols.gcid.create_csindex()
    uniqids_sorted = uniqids_out.copy(group, samplename, sortby='gcid', propindexes=True)

    h5.remove_node(group, samplename + '_working')
    h5.close()

    return grpstats


def write_dupstats(output, dupstats):
    w = csv.writer(output)
    w.writerow(('duplicated reads', 'groups of duplicated reads'))
    w.writerows(sorted(dupstats.items(), key=lambda x: x[0]))


def create_sqi_tables(h5, group, samplename, ncycles, nchannels):

    class SeqQual(tables.IsDescription):
        gcid = tables.Int64Col(pos=0)
        istart = tables.Int16Col(pos=1)
        seq = tables.StringCol(ncycles, pos=2)
        qual = tables.UInt8Col(shape=(ncycles,), pos=3)
        clones = tables.Int64Col(pos=4)

    subgroup = h5.create_group(group, samplename,
                               'SQI arrays for sample {}'.format(samplename))

    seqqual = h5.create_table(subgroup, 'seqqual', SeqQual,
                              'Sequences and quality scores for sample {}'.format(samplename))
    intensities = h5.create_earray(subgroup, 'intensities', tables.Int16Atom(),
                                   (0, ncycles, nchannels),
                                   'Signal intensities for sample {}'.format(samplename))

    return seqqual, intensities


def copy_unique_spots_for_tile(output_seqquals, output_intensities, seqfilters,
                               input_seqquals, input_intensities, iter_uniqids):
    iter_sqi = zip(input_seqquals.iterrows(), input_intensities.iterrows())

    if iter_uniqids is not None:
        uniqid_next = next(iter_uniqids)

    sqi_next = next(iter_sqi)
    uniq_passed = uniq_removed = total_written = 0

    while True:
        if sqi_next is None:
            sqi_next = next(iter_sqi, None)
            if sqi_next is None:
                break

        if iter_uniqids is not None:
            if uniqid_next is None:
                uniqid_next = next(iter_uniqids, None)

            if uniqid_next is None or sqi_next[0]['gcid'] != uniqid_next['gcid']:
                # a read removed during the deduplication
                sqi_next = None
                uniq_removed += 1
                continue

        uniq_passed += 1

        filterresults = [seqfilt(sqi_next[0]['seq']) for seqfilt in seqfilters]
        if all(filterresults):
            duplicity = uniqid_next['count'] if iter_uniqids is not None else 1
            output_seqquals.append([sqi_next[0].fetch_all_fields().tolist() + (duplicity,)])
            output_intensities.append([sqi_next[1]])
            total_written += 1

        sqi_next = uniqid_next = None

    if iter_uniqids is not None:
        stats = [('uniq', uniq_passed, uniq_removed)]
    else:
        stats = []
    stats.extend([seqfilt.get_stats() for seqfilt in seqfilters])
    stats.append(('all', total_written, uniq_passed + uniq_removed - total_written))

    return stats


def write_unique_spots(output_fn, sqiinput, uniqid_fn, seqfilters, samplename, tileid):
    if uniqid_fn is not None:
        uniqid_input = tables.open_file(uniqid_fn, mode='r')
    else:
        uniqid_input = None

    with tables.open_file(output_fn, mode='w') as output, \
         tables.open_file(sqiinput, mode='r') as inptbl:

        group = output.create_group('/', 'sqi', 'Deduplicated seq-qual-intensity arrays')

        seqquals = inptbl.get_node('/primary_sqi/{}/seqqual'.format(samplename))
        intensities = inptbl.get_node('/primary_sqi/{}/intensities'.format(samplename))

        if uniqid_input is not None:
            uniqids_table = uniqid_input.get_node('/uniqids/' + samplename)
            iter_uniqids = uniqids_table.where(tile_range_cond(tileid))
        else:
            iter_uniqids = None

        output_seqquals, output_intensities = (
                create_sqi_tables(output, group, samplename,
                                  intensities.shape[1], intensities.shape[2]))

        stats = copy_unique_spots_for_tile(output_seqquals, output_intensities,
                                           seqfilters, seqquals, intensities, iter_uniqids)

    if uniqid_input is not None:
        uniqid_input.close()

    return (tileid, stats)


def exec_prepare_fingerprints(options, infiles, fpregions, tmpdir):
    id_slices = [slice(begin-1, end) for begin, end in fpregions]

    with futures.ProcessPoolExecutor(options.parallel) as executor:
        jobs = []

        for tileid, infile in infiles.items():
            job = executor.submit(generate_id_sequences, infile, options.sample_name,
                                  tileid, id_slices, tmpdir)
            jobs.append(job)

        result_files = [job.result() for job in jobs]

    return result_files


def exec_sort_fingerprints_reduce_duplicates(options, fingerprint_files, tmpdir):
    # Sort sequence fingerprints and pick ones per duplicates.
    debug_output = open(options.debug, 'w') if options.debug is not None else None
    clonelist_out = open(options.clonelist, 'w') if options.clonelist is not None else None

    uniqid_filename = (options.output_uniqids if options.output_uniqids is not None
                       else os.path.join(str(tmpdir), 'uniqids.h5'))
    dupstats = sort_and_uniquify_duplicates(
            uniqid_filename, options.sample_name, fingerprint_files, options.parallel,
            debug_output, clonelist_out)

    if options.dupstats is not None:
        write_dupstats(open(options.dupstats, 'w'), dupstats)

    return uniqid_filename


def exec_filter_and_write_uniq_sequences(options, uniqid_filename, infiles, tmpdir):
    # Initialize sequence filters
    seqfilters = []
    for filterclass in (BalanceFilter, FixedUMIFilter):
        try:
            filterinst = filterclass(options)
        except:
            pass
        else:
            seqfilters.append(filterinst)

    # Create new SQI tables which only retain unique spots.
    temporary_sqi_table_name = os.path.join(tmpdir, '{tile}.h5')
    with futures.ProcessPoolExecutor(options.parallel) as executor:
        jobs = []

        for tileid, infile in infiles.items():
            job = executor.submit(write_unique_spots,
                                  temporary_sqi_table_name.format(tile=tileid),
                                  infile, uniqid_filename, seqfilters,
                                  options.sample_name, tileid)
            jobs.append(job)

        sqiwriting_counts = [job.result() for job in jobs]


    flatstats = []
    for tileid, tilestats in sqiwriting_counts:
        for filterstats in tilestats:
            flatstats.append((tileid,) + filterstats)

    stats = pd.DataFrame(flatstats, columns=['tile', 'filter', 'out', 'failed'])
    stats.insert(2, 'in', stats['out'] + stats['failed'])
    del stats['failed']

    return temporary_sqi_table_name, stats


def exec_merge_sqi_tables(output_fn, input_fn, tileids, samplename, totallength):
    filters = tables.Filters(complib='blosc:zlib')

    with tables.open_file(output_fn, 'w', filters=filters) as outtblf:
        sqigrp = outtblf.create_group('/', 'sqi')
        samplegrp = outtblf.create_group(sqigrp, samplename)

        seqquals_out = intensities_out = None
        written = 0

        for tileid in tileids:
            with tables.open_file(input_fn.format(tile=tileid), 'r') as tblf:
                seqquals = tblf.get_node('/sqi/{}/seqqual'.format(samplename))
                intensities = tblf.get_node('/sqi/{}/intensities'.format(samplename))

                if seqquals_out is None:
                    seqquals_out = outtblf.create_table(samplegrp, 'seqqual',
                                        seqquals.description, expectedrows=totallength)
                    intensities_out = outtblf.create_carray(samplegrp, 'intensities',
                                        tables.Int16Atom(),
                                        (totallength,) + intensities.shape[1:])

                writing_wnd_end = written + len(seqquals)

                seqquals_out.append(seqquals[:])
                intensities_out[written:writing_wnd_end] = intensities[:]

                written = writing_wnd_end

        if written != totallength:
            raise ValueError("Not all clusters have written in merging.")


def run(options, fpregions, infiles):
    with TemporaryDirectory() as tmpdir:
        if fpregions:
            fingerprint_files = exec_prepare_fingerprints(options, infiles, fpregions, tmpdir)
            uniqid_filename = exec_sort_fingerprints_reduce_duplicates(options,
                                    fingerprint_files, tmpdir)
        else:
            uniqid_filename = None

        # Apply sequence balance filters and write out filter-passed and unique sequences.
        temporary_sqi_table_name, filtstats = \
            exec_filter_and_write_uniq_sequences(options, uniqid_filename, infiles, tmpdir)

        totalwritten = filtstats[filtstats['filter'] == 'all']['out'].sum()

        # Merge it into one SQI table
        exec_merge_sqi_tables(options.output_uniq_sqi, temporary_sqi_table_name,
                              sorted(infiles.keys()), options.sample_name,
                              totalwritten)

        # Write the sequence filtration statistics
        totalcnt = filtstats.groupby('filter').agg('sum').reset_index()
        totalcnt['tile'] = '(total)'

        filtstats_full = pd.concat([
            totalcnt[['tile', 'filter', 'in', 'out']],
            filtstats.sort_values(by=['tile', 'filter']),
        ])
        filtstats_full['pass_rate'] = (
                filtstats_full['out'] / filtstats_full['in'] * 100.).round(3)
        filtstats_full.to_csv(options.filterstats, index=False)


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description=
                                     'Filters out duplicates and low diversity reads')
    parser.add_argument(dest='infiles', metavar='TILEID:FILE', type=str, nargs='+',
                        help='Paths to the SQI tables')
    parser.add_argument('--parallel', dest='parallel', metavar='N', type=int,
                        help='Number of parallel processors', default=4)
    parser.add_argument('--sample-name', dest='sample_name', metavar='NAME', type=str,
                        required=True, help='Number of the sample')

    parser.add_argument('--fingerprint-region', dest='regions', metavar='start:end', type=str,
                        default=[], action='append', help='Cycle range to be regarded')
    parser.add_argument('--umi-position', dest='umipos', metavar='N', type=int,
                        default=None, help='1-based position of UMI sequence')
    parser.add_argument('--umi-sequence', dest='umiseq', metavar='SEQ', type=str,
                        default=None, help='Sequence string to check for UMI')
    parser.add_argument('--umi-mismatch', dest='umimismatch', metavar='N', type=int,
                        default=1, help='Maximum mismatches to filter by UMI sequence')
    parser.add_argument('--balance-region', dest='balanceregion', metavar='start:end',
                        type=str, default=None,
                        help='1-based coordinate for interval for balance check')
    parser.add_argument('--balance-minimum', dest='balancemin', metavar='N', type=int,
                        default=None, help='Minimum count to filter in balance region')

    parser.add_argument('--output-uniqids', dest='output_uniqids', metavar='FILE', type=str,
                        default=None, help='Path to file for deduplicated list of cluster IDs')
    parser.add_argument('--output-uniq-sqi', dest='output_uniq_sqi', metavar='FILE', type=str,
                        required=True,
                        help='Path to file for deduplicated list of cluster IDs')
    parser.add_argument('--output-duplicity-stats', dest='dupstats', metavar='FILE', type=str,
                        help='Path to file for duplication statistics (CSV)', default=None)
    parser.add_argument('--output-filter-stats', dest='filterstats', metavar='FILE', type=str,
                        help='Path to file for filtering statistics (CSV)', default=None)
    parser.add_argument('--output-clone-list', dest='clonelist', metavar='FILE', type=str,
                        default=None, help='Path to file for clones list of duplicated reads')
    parser.add_argument('--output-debug', dest='debug', metavar='FILE', type=str,
                        default=None, help='Path to file for debug output')

    options = parser.parse_args()

    fpregions = [map(int, reg.split(':')) for reg in options.regions]
    sqiinputs = {int(spec.split(':', 1)[0]): spec.split(':', 1)[1]
                 for spec in options.infiles}

    return options, fpregions, sqiinputs


if __name__ == '__main__':
    options, fpregions, sqiinputs = parse_arguments()
    run(options, fpregions, sqiinputs)

