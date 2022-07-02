import numpy as np
import pytest
from epiread_tools.naming_conventions import *
from epiread_tools.em_utils import GenomicInterval
from epiread_tools.epireadToBedgraph import EpiToBedgraph
from epiread_tools.epiparser import EpireadReader, CoordsEpiread, CommentEpiSNP
import os

@pytest.fixture
def config_dict():
    epiread_config = {
      "genomic_intervals": ["chr1:100-200", "chr2:300-350", "chr2:400-650"],
      "cpg_coordinates": "tests/data/sample_cpg_file.bed.gz",
      "epiread_files": ["tests/data/fixed_withA.epiread.gz"],
      "outfile": "test_output.bedgraph",
      "epiformat": "old_epiread_A",
      "header": False,
      "bedfile": False,

    }
    return epiread_config


def test_bedgraph_from_intervals(config_dict):
    test_specific = {"genomic_intervals":["chr1:205499880-205500150","chr10:17599401-17600087"],
    "epiread_files" : ["tests/data/old_epiread_A_snps_with_comments.epiread.gz"]}
    config_dict.update(test_specific)
    runner = EpiToBedgraph(config_dict)
    runner.read_mixture()
    assert runner.matrices[0].shape == (15,3)

def test_snps(config_dict):
    test_specific = {"genomic_intervals": ["chr1:205499880-205500150"],
                     "epiread_files": ["tests/data/old_epiread_A_snps_with_comments.epiread.gz"]}
    config_dict.update(test_specific)
    parser = CommentEpiSNP(config_dict)
    mat, mapper = parser.file_list_to_csr("chr1", parser.genomic_intervals)
    assert mat.shape[0] == 15
    assert (mat[:,2141].toarray().flatten() == np.array([Y,A,C,C,A,A,Y,A,C,A] + [NOVAL]*5)).all()

def test_epiread_formats(config_dict):
    test_specific = {"genomic_intervals" : ["chr1:205511000-205511700"],
                     "epiread_files":  ["tests/data/problem_with_A.epiread.gz"]}
    config_dict.update(test_specific)
    parser_A = CoordsEpiread(config_dict)
    mat_A, mapper_A = parser_A.file_list_to_csr("chr1", parser_A.genomic_intervals)

    without_A = {"genomic_intervals" : ["chr1:205511000-205511700"],
                     "epiread_files":  ["tests/data/problem_without_A.epiread.gz"]}

    config_dict.update(without_A)
    parser = EpireadReader(config_dict)
    mat, mapper = parser.file_list_to_csr("chr1", parser.genomic_intervals)

    assert mat_A[:,mapper_A.abs_to_ind(205511553)].sum() == mat[:,mapper.abs_to_ind(205511553)].sum()

def test_overlapping_regions(config_dict):
    test_specific = {"genomic_intervals" : ["chr1:205511552-205511555", "chr1:205511552-205511555"],
    "epiread_files" :["tests/data/fixed_withA.epiread.gz"]}
    config_dict.update(test_specific)
    runner = EpiToBedgraph(config_dict)
    runner.read_mixture()
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

