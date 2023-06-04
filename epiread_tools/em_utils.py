###################################################
#
# Script: em_utils.py
# Author: Irene Unterman
# Description: useful classes and functions to
# handle WGBS data
#
# MIT License
#
# Copyright (c) 2022 irene unterman and ben berman
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ###################################################

from  epiread_tools.naming_conventions import *
import numpy as np
import pandas as pd
from itertools import chain
from bisect import bisect_right, bisect_left
import pysam
import re
from collections import defaultdict

class GenomicInterval:
    '''
    An obect of this class represents
     piece of the genome
     initialized with format chrN:start-stop
    '''

    def __init__(self, interval="chrN:0-0"):
        self.set_interval(interval)

    def set_interval(self, interval):
        '''
        set new interval
        :param interval: chrN:0-0 format
        :return: the GenomicInterval object
        '''
        self.interval = interval
        self.chrom, self.start, self.end = interval.replace(":", TAB).replace("-", TAB).split(TAB)
        self.start, self.end = int(self.start), int(self.end)
        return self

    def set_from_positions(self, chrom, start, end):
        '''
        set interval from chromosome, start and end
        :param chrom: chrN
        :param start: int
        :param end: int
        :return: the GenomicInterval object
        '''
        self.chrom, self.start, self.end = chrom, int(start), int(end)
        return self

    def __eq__(self, other):
        if not isinstance(other, GenomicInterval):
            # don't attempt to compare against unrelated types
            return NotImplemented

        return self.chrom == other.chrom and self.start == other.start and self.end == other.end

    def get_genomic_size(self):
        '''
        find genomic length of interval
        :return: int length
        '''
        return self.end - self.start

    def shift(self, n, inplace=False):
        '''
        shifts intervals by n bases, in place
        this does not take into account chromosome edges
        :param n: bases to move, neg in upstream
        :return:
        '''
        if inplace:
            self.start += n
            self.end += n
            return self
        else:
            return GenomicInterval().set_from_positions(self.chrom, self.start+n, self.end+n)

    def slop(self, n, inplace=False):
        '''
        this does not take into account chromosome edges
        :param n:  bases to add
        :param inplace: edit existing instance
        :return: slopped GenomicInterval
        '''
        start = max(self.start-n,0) #no neg coordinates
        end = self.end + n #might be bigger than chrom size, does it matter?
        if inplace:
            self.start = start
            self.end = end
            return self
        else:
            return GenomicInterval().set_from_positions(self.chrom, start, end)

    def __repr__(self):
        '''
        print interval in  chrN:0-0 format
        :return: str formatting
        '''
        return self.chrom+":"+str(self.start)+"-"+str(self.end)

    def __len__(self):
        return self.get_genomic_size()

def cut(fp, intervals):
    '''
    Tabix out interval from file, return reads
    in original format
    :param intervals: list of GenomicInterval of chunk
    :return: iterator of relevant records
    '''
    tabixfile = pysam.TabixFile(fp)
    for interval in intervals:
        try:
            yield from tabixfile.fetch(interval.chrom, interval.start, interval.end)
        except ValueError:  # no coverage for this region in file
            continue

class Mapper:
    '''
    keeps mapping from genomic to relative and vice versa
    '''

    def __init__(self,chrom, intervals, epiread_files, CpG_file, min_dist=0, slop=800):
        '''

        :param chrom: chromosome
        :param intervals: GenomicIntervals list
        :param epiread_files: list of epiread files
        :param CpG_file: file in bed format of cpg coordinates
        :param min_dist: minimun distance between intervals
        '''
        self.chrom = chrom
        self.original_intervals = intervals
        self.epiread_files = epiread_files
        self.CpG_file = CpG_file
        self.min_dist = min_dist
        self.slop = slop

        self.slopped_intervals = [x.slop(slop) for x in self.original_intervals]
        self.merged_intervals = self.merge_intervals(self.slopped_intervals)
        self.init_sample_ids()
        self.load_CpGs()
        self.load_snps()
        self.rel_intervals = self.get_nearest_rel_intervals(self.merged_intervals)
        self.init_rel_to_mat_ind()


    def load_snps(self, padding=2000):
        '''
        :param padding: how many bases to include before first
        interval and after last
        :return: mapping back and forth for absolute and
        relative snps
        '''
        #add padding to every interval, merge if adjacent, translate to ind
        bins = list(chain.from_iterable([(interval.start, interval.end) for interval in self.merged_intervals]))
        bins = np.sort(np.array(bins))
        bins[::2] -= padding
        bins[1::2] += padding
        bins = np.array(list(chain.from_iterable(merge_win_list(list(zip(bins[::2], bins[1::2]))))))#TODO: change, hard to read
        sizes = np.array([y-x for x, y in zip(bins[::2], bins[1::2])])
        cum_sizes = np.append(np.zeros(1), np.cumsum(sizes)).astype(int)
        max_snps = cum_sizes[-1]
        def abs_to_rel(abs_snp):
            index = np.digitize(abs_snp, bins) -1
            return cum_sizes[index//2] + abs_snp - bins[index]
        def rel_to_abs(rel_snp):
            index = np.digitize(rel_snp, cum_sizes) -1
            return bins[2*index] + rel_snp - cum_sizes[index]
        self.max_snps, self.snp_abs_to_rel, self.snp_rel_to_abs = max_snps, abs_to_rel, rel_to_abs

    def init_index_to_source(self, sources):
        '''
        initialize mapping of read indices
        to source sample
        :param sources: start and end read indices
        for each sample id (start, end, id)
        :return: func mapping read index to sample_id
        '''
        boundaries = [s[1] for s in sources]
        values = [s[2] for s in sources]
        def map_to_source(ind):
            return values[bisect_right(boundaries, ind)]
        self.index_to_source = map_to_source

    def init_sample_ids(self):
        '''
        initialize sample ID mapping
        :param epiread_files: list of samples
        :return: func mapping sample to ID
        '''
        sample_to_id = dict(zip(self.epiread_files,range(1,len(self.epiread_files)+1)))
        sample_to_id["All"] = 0
        self.sample_to_id = sample_to_id

    def get_sample_id_list(self):
        '''
        get list of all sample IDs
        :return: a list of sample ids
        '''
        return sorted(list(self.sample_to_id.values()))

    def load_CpGs(self):
        '''
        we want to load not only cpgs in intervals, but slop
        +- read length
        :return: mapping back and forth for absolute and
        relative coordinates
        '''
        CpGs = []
        for record in cut(self.CpG_file, self.merged_intervals):
            CpGs.append(int(record.split(TAB)[1]))  # start position in bed format
        CpGs = np.array(CpGs)
        rel_to_abs = dict(enumerate(CpGs))
        abs_to_rel = {v: k for k, v in rel_to_abs.items()}
        self.cpgs, self.abs_to_rel, self.rel_to_abs = CpGs, abs_to_rel, rel_to_abs

    def get_nearest_rel_intervals(self, genomic_intervals):
        '''
        from genomic coordinate, find nearest CpG
        :param genomic_intervals: genomic coordinates
        :return: relative coordinates
        '''
        rel_intervals = []
        for interval in genomic_intervals:
            rel_intervals.append((bisect_left(self.cpgs, interval.start), bisect_right(self.cpgs, interval.end)))
        return rel_intervals

    def init_rel_to_mat_ind(self):
        rel_cpgs = list(chain(*[list(range(x,y)) for (x,y) in self.rel_intervals]))
        index_to_relative = dict(enumerate(rel_cpgs))
        relative_to_index = {v: k for k, v in index_to_relative.items()}
        max_cpgs = len(rel_cpgs)
        def rel_to_ind(rel):
            return relative_to_index[rel]
        def ind_to_rel(ind):
            return index_to_relative[ind]
        def abs_to_ind(abs):
            return rel_to_ind(self.abs_to_rel[abs])
        def ind_to_abs(ind):
            return self.rel_to_abs[ind_to_rel(ind)]
        self.max_cpgs, self.abs_to_ind, self.ind_to_abs, self.rel_to_ind, self.ind_to_rel, self.all_rel = \
            max_cpgs, abs_to_ind, ind_to_abs, rel_to_ind, ind_to_rel, relative_to_index

    def get_ind_intervals(self, genomic_intervals):
        '''
        translate genomic coordinates to indices
        :param genomic_intervals: genomic coordinates
        :return: indices
        '''
        ind_intervals = []
        for rel_start, rel_end in self.get_nearest_rel_intervals(genomic_intervals):
            ind_intervals.append((self.rel_to_ind(rel_start), self.rel_to_ind(rel_end -1) + 1))
        return ind_intervals

    def merge_intervals(self, intervals):
        '''
        merge intervals if they are less than min dist apart
        :return: merged intervals
        '''
        #easier to load less chunks
        win_list = [(interval.start, interval.end) for interval in intervals]
        merged_list = merge_win_list(win_list, self.min_dist)
        new_intervals = [GenomicInterval().set_from_positions(self.chrom, x, y) for x, y in merged_list]
        return new_intervals




def merge_win_list(win_list, min_overlap=0):
    '''
    merge windows if they overlap
    :param win_list: list of windows
    :param min_overlap: minimal overlap for merging
    :return: merged windows
    '''
    new_windows = []
    current_start, current_end = win_list[0] #open first window
    for i, window in enumerate(win_list[1:]):
        win_start, win_end = window
        if win_start + min_overlap < current_end:
            current_end = max(current_end, win_end)
        else:
            new_windows.append((current_start, current_end))
            current_start, current_end = win_start, win_end
    new_windows.append((current_start, current_end))
    return new_windows


def split_intervals_to_chromosomes(genomic_intervals):
    '''
    find all chromosomes in intervals
    :return:
    '''
    by_chrom = defaultdict(list)
    for interval in genomic_intervals:
        by_chrom[interval.chrom].append(interval)
    return by_chrom

def bedgraph_to_intervals(fp, header=False):
    '''
    read bedgraph file and convert to GenomicInterval
    :param fp: file path to bedgraph, first cols are chrom, start, end
    :param header: does bedgraph file have header
    :return: list of GenomicInterval objects
    '''
    if header:
        df = pd.read_csv(fp, usecols=[0,1,2], skiprows=1, names=["chrom", "start", "end"], sep="\t")
    else:
        df = pd.read_csv(fp, usecols=[0,1,2], names=["chrom", "start", "end"], sep="\t")
    intervals = [GenomicInterval().set_from_positions(chrom, start, end) for chrom, start, end in df.to_records(index=False)]
    return intervals


def find_intersection(intervals, range_start, range_end):
    '''
    find intersection between intervals and range
    included book-ended eg: (3,3)
    :param intervals: np array of interval start, interval end
    :param range_start: range start
    :param range_end: range end
    :return: list of intersections
    '''
    intersection = []
    start_ind = np.searchsorted(intervals, range_start)
    end_ind = np.searchsorted(intervals, range_end, "right")
    if start_ind %2: #even:
       start_ind -=1
    if end_ind%2: #odd
        end_ind += 1
    for interval_start, interval_end in zip(intervals[start_ind:end_ind:2], intervals[start_ind+1:end_ind:2]):
        intersection.append((max(range_start, interval_start), min(range_end, interval_end)))
    return intersection

def calc_coverage(matrices):
    '''
    calculate coverage per CpG
    :return:np array of coverage sums
    '''
    coverage = [np.squeeze(np.asarray((methylation_matrix == METHYLATED).sum(axis=0)))+ \
                    np.squeeze(np.asarray((methylation_matrix == UNMETHYLATED).sum(axis=0)))
                     for methylation_matrix in matrices]
    coverage = [x if type(x)==np.ndarray else np.array([x]) for x in coverage]
    return coverage

def calc_methylated(matrices):
    '''
    calc sum methylated per item in list
    :param matrices: list of np arrays
    :return: np array of methylation sums
    '''
    methylation = [np.squeeze(np.asarray((methylation_matrix == METHYLATED).sum(axis=0)))
                        for methylation_matrix in matrices]
    return methylation

def calc_percent_U(mat, u_threshold=0.25):
    '''
    calc percent unmehtylated reads in mat
    DOES NOT FILTER BY #CPG
    DOES NOT FILTER EMPTY READS
    :param mat: np array
    :param u_threshold: maximal prop of methylation
    :return: %U reads
    '''
    meth = (mat == METHYLATED).sum(axis=1)
    unmeth = (mat == UNMETHYLATED).sum(axis=1)
    is_u = (meth / (unmeth + meth)) <= u_threshold
    percent_u = np.sum(is_u) / mat.shape[0]
    return percent_u

def in_intervals(cpg, intervals):
    '''
    check if cpg in intervals. assumes cpg and
    all intervals are on the same chromosome
    :param cpg: start position, e.g. 21244961
    :param intervals: list of Genomic intervals
    :return: cpg in intervals
    '''
    sorted_intervals = sorted(intervals, key=lambda x: x.start)
    if cpg < sorted_intervals[0].start:
        return False
    for interval in sorted_intervals:
        if cpg > interval.start and cpg < interval.end:
            return True
    return False

def jsonconverter(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, datetime.datetime):
        return obj.__str__()