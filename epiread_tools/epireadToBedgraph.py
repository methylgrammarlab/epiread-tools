###################################################
#
# Script: epireadToBedGraph.py
# Author: Irene Unterman and Ben Berman
# Description: parse epiread files and save mean
# methylation per position
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
###################################################
import json

import click
import sys

import numpy as np
from epiread_tools.epiparser import EpireadReader, CoordsEpiread, PatReader
from epiread_tools.epiformat import epiformat_to_reader
from epiread_tools.naming_conventions import *
from epiread_tools.em_utils import calc_coverage, calc_methylated

class EpiToBedgraph:

    def __init__(self, config):
        '''
        Initializes an instance of epiread to bedgraph transformer

        Args:
            config (dict): Configuration settings for the EpiToBedgraph instance.
            Should contain:
            epiformat - epiread format (old_epiread, old_epiread_A, pat)
            outfile - output file path
        '''
        self.config = config
        self.reader = epiformat_to_reader[self.config["epiformat"]]

    def read_mixture(self):
        '''
        Reads the epiread and retrieves matrices, CpGs, and origins.
        '''
        reader = self.reader(self.config)
        self.interval_order, self.matrices, self.cpgs, self.origins = reader.get_matrices_for_intervals()
        print(f"Read {len(self.interval_order)} intervals from {self.config['epiread_files']}")
#        print(f"\t{self.interval_order}\n")
        print(f"Read {len(self.matrices)} matrices from {self.config['epiread_files']}\n")
#        print(f"\t{self.matrices}\n")
        print(f"Read {len(self.cpgs)} CpGs from {self.config['cpg_coordinates']}\n")
#        print(f"\t{self.cpgs}\n")
        print(f"Read {len(self.origins)} origins from {self.config['epiread_files']}\n")
#        print(f"\t{self.origins}\n")


    def calc_coverage(self):
        '''
        calculate coverage per CpG
        '''
        self.coverage = calc_coverage(self.matrices)

    def calc_methylated(self):
        '''
        Calculates the methylated counts per CpG based on the matrices.
        '''
        self.methylation = calc_methylated(self.matrices)

    def calc_mean_methylation(self):
        '''
        mean methylation per CpG with coverage
        '''
        #skip all-zero columns, save indices for later
        self.mean_methylation = []
        for i in range(len(self.matrices)):
            mean_methylation = np.divide(self.methylation[i],self.coverage[i],where=self.coverage[i]>0)
            self.mean_methylation.append(mean_methylation)

    def write_bedgraph(self):
        '''
        append output to output file
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
                print(f"skipped {len(self.coverage[i])-len(arr)} zero coverage CpGs")
            elif len(self.coverage[i]) == 1 and not self.coverage[i][0]:
                arr = np.array([]) #no coverage
                print(f"skipped 1 zero coverage CpG")

            if arr.any():
                out_arrs.append(arr)

        if out_arrs:
            with open(self.config["outfile"], "w") as outfile:
                np.savetxt(outfile, np.vstack(out_arrs), delimiter="\t", fmt="%s")
        else:
            print("No data to write to the output file.")

    def tobedgraph(self):
        self.read_mixture()
        self.calc_coverage()
        self.calc_methylated()
        self.calc_mean_methylation()
        self.write_bedgraph()

#%%


@click.command(context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.option('--cpg_coordinates', help='sorted cpg bed file')
@click.option('--epiread_files', help='comma delimited epiread paths (required)')
@click.option('--outfile', help='output file path')
@click.option('-j', '--json', help='run from json config file')
@click.option('-i', '--genomic_intervals', help='interval(s) to process. formatted chrN:start-end, separated by commas')
@click.option('-b', '--bedfile', help='bed file chrom start end with interval(s) to process. tab delimited',
              is_flag=True, default=False)
@click.option('--header', is_flag=True, default=False, help="bedgraph with regions to process has header")
@click.option('-A', '--coords', is_flag=True, help='epiread files contain coords', default=False)
@click.option('--epiformat',
              type=click.Choice(['old_epiread','old_epiread_A','pat'], case_sensitive=False))
@click.version_option()
@click.pass_context
def main(ctx, **kwargs):
    """ biscuit epiread to bedgraph converter. any command line options will override config"""

    # I merged in new parsing code, but haven't tested it yet.
    # Commented out below is the old code that worked.

        # if kwargs["json"]:
        #     with open(kwargs["json"], "r") as json_file:
        #         config = json.load(json_file)
        # else:
        #     config= {"epiformat":"old_epiread", "bedfile":False}

        # for key, value in kwargs.items():
        #     config.setdefault(key, value)

        # config.update(dict([item.strip('--').split('=') for item in ctx.args]))
        # if config["epiread_files"]:
        #     if not isinstance(config["epiread_files"], list):
        #         config['epiread_files'] = config['epiread_files'].split(",")
        # else:
        #     print("No epiread files provided.\n\n")
        #     print(ctx.get_help())  # Print command-line usage        
        #     sys.exit(1)

        # if config["intervals"]:
        #     config['genomic_intervals'] = config["intervals"].split(",")
        # if config["coords"]: #flag
        #     config['epiformat'] = "old_epiread_A"


    config = {"epiformat":"old_epiread"}
    config.update(kwargs)
    config.update(dict([item.strip('--').split('=') for item in ctx.args]))
    if kwargs["json"] is not None:
        with open(kwargs["json"], "r") as jconfig:
            config.update(json.load(jconfig))



    if config["epiread_files"]:
        config['epiread_files'] = config['epiread_files'].split(",")
    if not config["bedfile"]:
        config['genomic_intervals'] = config["intervals"].split(",")
    if config["coords"]: #flag
        config['epiformat'] = "old_epiread_A"
    runner = EpiToBedgraph(config)
    runner.tobedgraph()


if __name__ == '__main__':
    main()
# config = {"cpg_coordinates": "/Users/ireneu/berman_lab/ALS/hg19_pat_cpg.bed.gz", "bedfile":True,
#           "genomic_intervals":"/Users/ireneu/berman_lab/ALS/pat_U250.bed", "outfile":"/Users/ireneu/berman_lab/ALS/test.bedgraph",
#           "epiformat":"pat", "header":False, "epiread_files":["/Users/ireneu/berman_lab/ALS/Small-Int-Ep_atlas_epipaths.epiread.gz"]}
# runner = EpiToBedgraph(config)
# runner.tobedgraph()