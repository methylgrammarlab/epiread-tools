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


class BedgraphRunner():

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


    def parse_reads(self, chrom, intervals):
        '''
        get data amd mapper
        :return:
        '''
        parser = Parser(chrom, intervals, self.epiread_files, self.cpg_locations, self.epiformat)
        self.methylation_matrix, self.mapper = parser.parse()


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
        with open(self.outfile, "a+") as outfile:
            np.savetxt(outfile, output_array, delimiter=TAB, fmt='%s')

    def tobedgraph(self):
        for chrom, intervals in self.intervals_per_chrom.items():
            self.parse_reads(chrom, intervals)
            self.calc_coverage()
            self.calc_mean_methylation()
            self.write_bedgraph(chrom)

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
    """ input cpg_coordinates file, epiread files separated by comma and output file"""
    click.echo("biscuit epiread to bedgraph converter")
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
    runner = BedgraphRunner(genomic_intervals, cpg_coordinates, epiread_files, outfile, epiformat,
                            kwargs["header"], kwargs["bedfile"])
    runner.tobedgraph()


if __name__ == '__main__':
    main()
