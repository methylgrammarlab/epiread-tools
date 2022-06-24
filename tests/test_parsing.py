import pytest
from epiread_tools.naming_conventions import *
from epiread_tools.epireadToBedgraph import BedgraphRunner



def test_bedgraph_from_intervals(self):
    genomic_intervals=["chr1:205499880-205500511"]
    cpg_coordinates = "./tests/data/sample_cpg_file.bed.gz"
    epiread_files = ["./tests/data/old_epiread_A_snps_with_comments.epiread.gz"]
    runner = BedgraphRunner(genomic_intervals, cpg_coordinates, epiread_files, outfile=None, header=False, bedfile=False)
    for chrom, intervals in runner.intervals_per_chrom.items():
        runner.parse_reads(chrom, intervals)
    assert runner.methylation_matrix.shape == (15,3)

def test_bedgraph_from_bedfile(self):
    pass
    #with & without header

def test_old_epiread(self):
    pass

def test_coord_epiread(self):
    pass

def test_mapper(self):
    pass

