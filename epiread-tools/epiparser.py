###################################################
#
# Script: epiparser.py
# Author: Irene Unterman
# Description: parse epiread files for EM
###################################################

import numpy as np
import pysam
from naming_conventions import *
from em_utils import asEpiread, Mapper
import scipy.sparse as sp
from itertools import tee
#%%


def parse_epireads(chrom, intervals, epiread_files, CpG_file, one_based, parse_snps=True):
    '''
    read all epiread files to sparse matrix, optionally parse mutation data
    :param chrom: chromosome
    :param intervals: list of GenomicInterval of chunk
    :param epiread_files: list of input epireads
    :param CpG_file
    :param parse_snps: add snp data
    :return: matrix of methylation, matrix of SNPs (if parse_snps),
    mapping from abs to rel and back
    '''
    mapper = Mapper(chrom, intervals, epiread_files, CpG_file, one_based) #init mapping
    small_matrices = []
    if parse_snps:
        snp_matrices = []
    sources = []
    i = 0
    for epi_file in epiread_files:
        if parse_snps:
            epi_iter, snp_iter = tee(cut_epiread(intervals, epi_file))
        else:
            epi_iter = cut_epiread(intervals, epi_file)
        small_matrix = epiread_to_csr(epi_iter, mapper.rel_intervals, mapper.abs_to_rel,
                                      mapper.rel_to_ind, mapper.max_cpgs)
        small_matrices.append(small_matrix)
        sources.append((i, i + small_matrix.shape[0], mapper.sample_to_id[epi_file]))
        i += small_matrix.shape[0]
        #handle snps
        if parse_snps:
            snp_matrices.append(snp_to_csr(snp_iter, mapper.snp_abs_to_rel, small_matrix.shape[0],
                                           mapper.max_snps))
    mapper.init_index_to_source(sources)
    methylation_matrix, snp_matrix = sp.vstack(small_matrices, dtype=int), None
    if parse_snps:
        snp_matrix = sp.vstack(snp_matrices, dtype=int)
    return methylation_matrix, snp_matrix, mapper


def cut_epiread(intervals, epiread_file):
    '''
    Tabix out interval from file, return reads
    in original format
    :param intervals: list of GenomicInterval of chunk
    :param epiread_file: one Tabix-indexed file in epiread format
    :return: iterator of relevant records
    '''
    tabixfile = pysam.TabixFile(epiread_file)
    for interval in intervals:
        try:
            yield from tabixfile.fetch(interval.chrom, interval.start, interval.end)
        except ValueError: #no coverage for this region in file
            # with open("coverage.txt", "a+") as outfile:
            #     outfile.write(epiread_file+TAB+str(interval)+"\n")
            continue
#%%
def find_intersection(intervals, range_start, range_end):
    '''
    find intersection between intervals and range
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



def epiread_to_csr(epiread_iterator, rel_intervals, abs_to_rel, rel_to_ind, max_CpGs):
    '''
    map reads from single source to sparse matrix
    :param epiread_iterator: iterator of epiread records
    :param rel_intervals: relative CpG coordinates to parse
    :param abs_to_rel: dict from genomic coordinates to relative
    :param rel_to_ind: func mapping relative coordinates to index in matrix
    :param max_CpGs: total number of CpGs in matrix
    :return: sparse matrix of reads
    '''
    row = []
    col = []
    data = []
    n_reads = 0
    rel_intervals = np.array(rel_intervals).flatten()
    for i, epiread in enumerate(epiread_iterator):
        n_reads = i
        record = asEpiread(*epiread.split(TAB))
        rel_start = abs_to_rel[record.get_start()]
        rel_end = abs_to_rel[record.get_end()] #TODO: implement!
        for intersect_start, intersect_end in find_intersection(rel_intervals, rel_start, rel_end):
            if not intersect_end - intersect_start: #overlap length is 0
                continue
            row.extend([i]*(intersect_end-intersect_start))
            col.extend(list(range(rel_to_ind(intersect_start), rel_to_ind(intersect_end -1) + 1)))
            for cpg in record.methylation[intersect_start-rel_start:intersect_end-rel_start]:
                data.append(methylation_state[cpg])
    return sp.csr_matrix((data, (row, col)), shape=(n_reads+1, max_CpGs), dtype=int)



def snp_to_csr(epiread_iterator, snp_abs_to_rel, n_reads, max_snps):
    '''
    map snps from single source to sparse matrix
    :param epiread_iterator: iterator of epiread records
    :param snp_abs_to_rel: dict from genomic coordinates to relative
    :param n_reads: number of methylation reads in iterator
    :param max_snps: maximal number of snps
    :return: sparse matrix of reads
    '''
    row= []
    col = []
    data = []
    for i, epiread in enumerate(epiread_iterator):
        record = asEpiread(*epiread.split(TAB))
        if record.has_snps():
            rel_start = snp_abs_to_rel(record.get_snp_start())
            for dist, snp in record.get_snps():
                if rel_start+int(dist) < max_snps and rel_start+int(dist): #doesn't pass interval
                    row.append(i)
                    col.append(rel_start+int(dist))
                    data.append(dna_to_ind[snp])
    return sp.csr_matrix((data, (row, col)), shape=(n_reads, max_snps), dtype=int)



