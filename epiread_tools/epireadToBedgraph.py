###################################################
#
# Script: epireadToBedGraph.py
# Author: Irene Unterman
# Description: parse epiread files and save mean
# methylation per position
###################################################

import click
import sys

import numpy as np
from  epiread_tools.epiparser import Parser
from  epiread_tools.em_utils import GenomicInterval, split_intervals_to_chromosomes, bedgraph_to_intervals
from epiread_tools.naming_conventions import *


class EpiRunner():

    def __init__(self, genomic_intervals, cpg_locations, epiread_files, outfile, epiformat, header=False, bedfile=True):
        self.header = header
        if bedfile:
            self.genomic_intervals = bedgraph_to_intervals(genomic_intervals, self.header)
        else:
            self.genomic_intervals = [GenomicInterval(x) for x in genomic_intervals]

        self.original_intervals = genomic_intervals
        self.cpg_locations = cpg_locations
        self.epiread_files = epiread_files
        self.outfile = outfile
        self.epiformat = epiformat
        self.intervals_per_chrom = split_intervals_to_chromosomes(self.genomic_intervals)
        self.matrices = []
        self.interval_order = []
        self.cpgs = []


    def parse_one_chromosome(self, chrom, intervals):
        '''
        get data amd mapper
        :return:
        '''
        intervals = sorted(intervals, key=lambda x: x.start)
        self.interval_order.extend(intervals)
        parser = Parser(chrom, intervals, self.epiread_files, self.cpg_locations, self.epiformat)
        methylation_matrix, mapper = parser.parse()
        window_list = mapper.get_ind_intervals(intervals)
        for start, end in window_list:
            slice = methylation_matrix[:,start:end]
            self.matrices.append(slice[slice.getnnz(1)>0]) #remove empty rows
            self.cpgs.append(np.array([mapper.ind_to_abs(x) for x in range(start, end)]))


    def parse_multiple_chromosomes(self):
        for chrom, intervals in self.intervals_per_chrom.items():
            self.parse_one_chromosome(chrom, intervals)


    def calc_coverage(self):
        '''
        calculate coverage per CpG
        :return:
        '''
        self.coverage = [np.squeeze(np.asarray((methylation_matrix == METHYLATED).sum(axis=0)))+ \
                        np.squeeze(np.asarray((methylation_matrix == UNMETHYLATED).sum(axis=0)))
                         for methylation_matrix in self.matrices]

    def calc_methylated(self):
        self.methylation = [np.squeeze(np.asarray((methylation_matrix == METHYLATED).sum(axis=0)))
                            for methylation_matrix in self.matrices]

    def calc_mean_methylation(self):
        '''
        mean methylation per CpG with coverage
        :return:
        '''
        #skip all-zero columns, save indices for later
        self.mean_methylation = []
        for i in range(len(self.matrices)):
            mean_methylation = np.divide(self.methylation[i],self.coverage[i],where=self.coverage[i]>0)
            self.mean_methylation.append(mean_methylation)

    def write_bedgraph(self):
        '''
        append output to output file
        :return:
        '''
        out_arrs = []
        for i, interval in enumerate(self.interval_order):
            arr=np.zeros(shape=(len(self.cpgs[i]), 5), dtype=object)
            arr[:,0]=interval.chrom
            arr[:,1]=self.cpgs[i]
            arr[:,2]=self.cpgs[i] + 1
            arr[:,3]=self.mean_methylation[i]
            arr[:,4]=self.coverage[i]
            arr = arr[self.coverage[i]>0,:]
            out_arrs.append(arr)

        with open(self.outfile, "w") as outfile:
            np.savetxt(outfile, np.vstack(out_arrs), delimiter=TAB, fmt='%s')

    def tobedgraph(self):
        self.parse_multiple_chromosomes()
        self.calc_coverage()
        self.calc_methylated()
        self.calc_mean_methylation()
        self.write_bedgraph()

#%%

@click.command()
@click.argument('cpg_coordinates')
@click.argument('epireads')
@click.argument('outfile')
@click.option('-i', '--intervals', help='interval(s) to process. formatted chrN:start-end, separated by commas')
@click.option('-b', '--bedfile', help='bed file chrom start end with interval(s) to process. tab delimited',
              is_flag=True, default=False)
@click.option('--header', is_flag=True, help="bedgraph with regions to process has header")
@click.option('-A', '--coords', is_flag=True, help='epiread files contain coords', default=False)
def main(cpg_coordinates, epireads, outfile ,*args, **kwargs):
    """ biscuit epiread to bedgraph converter"""
    epiread_files = epireads.split(",")

    if kwargs["intervals"]:
        genomic_intervals = kwargs["intervals"].split(",")
    elif kwargs["bedfile"]:
        genomic_intervals = kwargs["bedfile"]
    else:
        raise ValueError("either specify intervals or add bed file. For whole genome use -b with chrom sizes")
    epiformat = "old_epiread"
    if kwargs["coords"]:
        epiformat = "old_epiread_A"
    runner = EpiRunner(genomic_intervals, cpg_coordinates, epiread_files, outfile, epiformat,
                       kwargs["header"], kwargs["bedfile"])
    runner.tobedgraph()


if __name__ == '__main__':
    main()
#%%

# genomic_intervals=["chr1:205511000-205511700"]
# cpg_coordinates = "/Users/ireneu/PycharmProjects/epiread-tools/tests/data/sample_cpg_file.bed.gz"
# with_A = ["/Users/ireneu/PycharmProjects/epiread-tools/tests/data/problem_with_A.epiread.gz"]
# without_A = ["/Users/ireneu/PycharmProjects/epiread-tools/tests/data/problem_without_A.epiread.gz"]
# new_runner = EpiRunner(genomic_intervals, cpg_coordinates, with_A, outfile=None, epiformat="old_epiread_A",
#                    header=False, bedfile=False)
# new_runner.parse_multiple_chromosomes()
# old_runner = EpiRunner(genomic_intervals, cpg_coordinates, without_A, outfile=None, epiformat="old_epiread",
#                    header=False, bedfile=False)
# old_runner.parse_multiple_chromosomes()