import numpy as np
import pytest
from epiread_tools.naming_conventions import *
from epiread_tools.em_utils import GenomicInterval
from epiread_tools.epireadToBedgraph import EpiRunner
from epiread_tools.epiparser import Parser
import os



def test_bedgraph_from_intervals():
    genomic_intervals=["chr1:205499880-205500150"]
    cpg_coordinates = "tests/data/sample_cpg_file.bed.gz"
    epiread_files = ["tests/data/old_epiread_A_snps_with_comments.epiread.gz"]
    runner = EpiRunner(genomic_intervals, cpg_coordinates, epiread_files, outfile=None, epiformat="old_epiread_A",
                       header=False, bedfile=False)
    runner.parse_multiple_chromosomes()
    assert runner.matrices[0].shape == (15,3)

def test_snps():
    chrom = "chr1"
    intervals = [GenomicInterval("chr1:205499880-205500150")]
    cpg_file = "tests/data/sample_cpg_file.bed.gz"
    epiread_files = ["tests/data/old_epiread_A_snps_with_comments.epiread.gz"]
    epi_format = "snps_with_comments"
    parser = Parser(chrom, intervals, epiread_files, cpg_file, epi_format)
    mat, mapper = parser.parse()
    assert mat.shape[0] == 15
    assert (mat[:,2141].toarray().flatten() == np.array([Y,A,C,C,A,A,Y,A,C,A] + [NOVAL]*5)).all()

def test_epiread_formats():
    intervals=[GenomicInterval("chr1:205511000-205511700")]
    cpg_coordinates = "tests/data/sample_cpg_file.bed.gz"
    with_A = ["tests/data/problem_with_A.epiread.gz"]
    without_A = ["tests/data/problem_without_A.epiread.gz"]
    parser_A = Parser("chr1", intervals, with_A, cpg_coordinates, "old_epiread_A")
    mat_A, mapper_A = parser_A.parse()
    parser = Parser("chr1", intervals, without_A, cpg_coordinates, "old_epiread")
    mat, mapper = parser.parse()

    assert mat_A[:,mapper_A.abs_to_ind(205511553)].sum() == mat[:,mapper.abs_to_ind(205511553)].sum()

def test_overlapping_regions():
    genomic_intervals = ["chr1:205511552-205511555", "chr1:205511552-205511555"]
    cpg_coordinates = "tests/data/sample_cpg_file.bed.gz"
    epiread = ["tests/data/fixed_withA.epiread.gz"]
    runner = EpiRunner(genomic_intervals, cpg_coordinates, epiread,
                           outfile=None,
                           epiformat="old_epiread_A",
                           header=False, bedfile=False)
    runner.parse_multiple_chromosomes()
    runner.calc_coverage()
    assert runner.coverage[0]==runner.coverage[1]
    assert runner.coverage[0]==16



def test_bedgraph_from_bedfile():
    pass
    #with & without header

def test_old_epiread():
    pass

def test_coord_epiread():
    pass

def test_mapper():
    pass

