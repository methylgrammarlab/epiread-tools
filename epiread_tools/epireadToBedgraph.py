###################################################
#
# Script: epireadToBedGraph.py
# Author: Irene Unterman
# Description: parse epiread files and save mean
# methylation per position
###################################################
import json

import click
import sys

import numpy as np
from  epiread_tools.epiparser import EpireadReader, CoordsEpiread, epiformat_to_reader
from epiread_tools.naming_conventions import
from epiread_tools.em_utils import calc_coverage, calc_methylated


class EpiToBedgraph():

    def __init__(self, config):
        self.config = config
        self.reader = epiformat_to_reader[self.config["epiformat"]]

    def read_mixture(self):
        reader = self.reader(self.config)
        self.interval_order, self.matrices, self.cpgs = reader.get_matrices_for_intervals()

    def calc_coverage(self):
        '''
        calculate coverage per CpG
        :return:
        '''
        self.coverage = calc_coverage(self.matrices)

    def calc_methylated(self):
        self.methylation = calc_methylated(self.matrices)

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
            if len(self.coverage[i]) > 1: #doesn't work when length == 1
                arr = arr[(self.coverage[i]>0),:]
            elif len(self.coverage[i]) == 1 and not self.coverage[i][0]:
                arr = np.array([]) #no coverage
            if arr.any():
                out_arrs.append(arr)

        with open(self.outfile, "w") as outfile:
            np.savetxt(outfile, np.vstack(out_arrs), delimiter=TAB, fmt='%s')

    def tobedgraph(self):
        self.read_mixture()
        self.calc_coverage()
        self.calc_methylated()
        self.calc_mean_methylation()
        self.write_bedgraph()

#%%

@click.command()
@click.option('--cpg_coordinates', help='sorted cpg bed file')
@click.option('--epireads', help='comma delimited epiread paths')
@click.option('--outfile', help='output file path')
@click.option('-j', '--json', help='run from json config file')
@click.option('-i', '--intervals', help='interval(s) to process. formatted chrN:start-end, separated by commas')
@click.option('-b', '--bedfile', help='bed file chrom start end with interval(s) to process. tab delimited',
              is_flag=True, default=False)
@click.option('--header', is_flag=True, help="bedgraph with regions to process has header")
@click.option('-A', '--coords', is_flag=True, help='epiread files contain coords', default=False)
@click.version_option()
def main(**kwargs):
    """ biscuit epiread to bedgraph converter. any command line options will override config"""
    if kwargs["json"]:
        config = json.load(kwargs["json"])
    else:
        config= {"epiformat":"old_epiread", "bedfile":False}
    config.update(kwargs)
    if kwargs['epireads']:
        config['epireads'] = kwargs['epireads'].split(",")
    if kwargs["intervals"]:
        config['genomic_intervals'] = kwargs["intervals"].split(",")
    if kwargs["coords"]:
        config['epiformat'] = "old_epiread_A"
    if not config['genomic_intervals'] and not config['bedfile']:
        raise ValueError("either specify intervals or add bed file. For whole genome use -b with chrom sizes")
    runner = EpiToBedgraph(config)
    runner.tobedgraph()


if __name__ == '__main__':
    main()

