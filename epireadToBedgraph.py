###################################################
#
# Script: epireadToBedGraph.py
# Author: Irene Unterman
# Description: parse epiread files and save mean
# methylation per position
###################################################

import sys
import argparse
import os
import numpy as np
from collections import defaultdict
from epiparser import parse_epireads
from em_utils import GenomicInterval, Mapper, split_intervals_to_chromosomes, bedgraph_to_intervals
from naming_conventions import *
import pandas as pd
import scipy.sparse as sp
import json

class BedgraphRunner():

    def __init__(self, genomic_intervals, cpg_locations, epiread_files, outfile, header=False, bedgraph=True):
        self.header = header
        if bedgraph:
            self.genomic_intervals = bedgraph_to_intervals(genomic_intervals, self.header)
        else:
            # self.genomic_intervals = [GenomicInterval().set_from_positions(chrom, start, end) for chrom, start, end in genomic_intervals]
            self.genomic_intervals = [GenomicInterval(x) for x in genomic_intervals]

        self.original_intervals = genomic_intervals
        self.cpg_locations = cpg_locations
        self.epiread_files = epiread_files
        self.outfile = outfile
        self.intervals_per_chrom = split_intervals_to_chromosomes(self.genomic_intervals)


    def parse_reads(self, chrom, intervals):
        '''
        set data amd mapper
        :return:
        '''
        self.methylation_matrix, _, self.mapper = parse_epireads(chrom, intervals,
                                    self.epiread_files, self.cpg_locations,False, False)

    def calc_coverage(self):
        '''
        calculate coverage per CpG
        :return:
        '''
        self.coverage = np.squeeze(np.asarray((self.methylation_matrix == METHYLATED).sum(axis=0)))+ \
                        np.squeeze(np.asarray((self.methylation_matrix == UNMETHYLATED).sum(axis=0)))

    def calc_methylated(self):
        self.methylation = np.squeeze(np.asarray((self.methylation_matrix == METHYLATED).sum(axis=0)))

    def calc_mean_methylation(self):
        '''
        mean methylation per CpG with coverage
        :return:
        '''
        methylated_per_col = np.squeeze(np.asarray((self.methylation_matrix == METHYLATED).sum(axis=0)))
        #skip all-zero columns, save indices for later
        self.mean_methylation = np.divide(methylated_per_col[self.coverage>0],self.coverage[self.coverage>0])

    def write_bedgraph(self, chrom):
        '''
        append chromosome output to output file
        :param chrom: chromosome name
        :return:
        '''
        col_indices = np.arange(len(self.coverage))[self.coverage>0]
        output_array = np.zeros(shape=(len(col_indices), 5), dtype=object)
        for i in range(len(col_indices)):
            output_array[i, :] = chrom, self.mapper.ind_to_abs(col_indices[i]), \
                                 self.mapper.ind_to_abs(col_indices[i])+ 1, self.mean_methylation[i], \
                                 self.coverage[col_indices[i]]
        with open(self.outfile, "a") as outfile:
            np.savetxt(outfile, output_array, delimiter=TAB, fmt='%s')

    def tobedgraph(self):
        #first split by chromosome
        self.intervals_per_chrom = split_intervals_to_chromosomes(self.genomic_intervals)
        for chrom, intervals in self.intervals_per_chrom.items():
            self.parse_reads(chrom, intervals)
            self.calc_coverage()
            self.calc_mean_methylation()
            self.write_bedgraph(chrom)
#%%
parser = argparse.ArgumentParser(description='Transform epiread file(s) to bedgraph')
parser.add_argument('cpg_coordinates', help='file with coordinates of CpGs')
parser.add_argument('epireads', help='file(s) to process. gziiped and indexed. multiple files will be aggregated, '+\
                                     'separated by commas')
parser.add_argument('outfile', help='path for output file')
parser.add_argument('-i','--intervals', help='interval(s) to process. formatted chrN:start-end, separated by commas',
                    default=False)
parser.add_argument('-b','--bedfile', help='bed file chrom start end with interval(s) to process. tab delimited',
                    default=False)
parser.add_argument('-g','--genome', help='file with chromosome:bp, used to parse entire file if interval not specified',
                    default=False)
parser.add_argument('-d', '--header', action='store_true', help="bedgraph with regions to process has header")

args = parser.parse_args()
epiread_files = args.epireads.split(",")

if args.intervals:
    genomic_intervals = args.intervals.split(",")
elif args.bedfile:
    genomic_intervals = args.bedfile
else:
    pass
    #TODO: implement genome option
runner = BedgraphRunner(genomic_intervals, args.cpg_coordinates, epiread_files, args.outfile, args.header, args.bedfile)
runner.tobedgraph()


#%%
# genomic_intervals=["chr1:205598000-205603169"]
# cpg_coordinates = "/Users/ireneu/PycharmProjects/dmr-cleanup/in-silico/epireads/hg19.CpG.bed.gz"
# epiread_files = ["/Users/ireneu/PycharmProjects/dmr-cleanup/epiread_format/small_Pancreas-Beta-Z0000043H.after_fix_bug_dash.epiread.gz"]
# outfile="."
# header=False
# bedfile=False
# runner = BedgraphRunner(genomic_intervals, cpg_coordinates, epiread_files, outfile, header, bedfile)
# runner.tobedgraph()