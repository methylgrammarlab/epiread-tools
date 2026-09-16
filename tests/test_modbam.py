###################################################
#
# Script: modbam_parser.py
# Author: Ben Berman
# Description: parse epiread files for EM
#
# Copyright (c) 2023 ben berman, base4analytics, and volitionRx
# ###################################################

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
      "cpg_coordinates": "tests/data/epiread_cpg_file.bed.gz",
      "epiread_files": ["tests/data/old_epiread.epiread.gz"],
      "outfile": "test_output.bedgraph",
      "epiformat": "old_epiread",
      "header": False,
      "bedfile": False,
    }
    return epiread_config

def test_modbam(config_dict):
    test_specific = {"genomic_intervals": ["chr16:616000-617000","chr16:19113000-19116000"],
    "epiread_files" : ["tests/data/test-fw-rev.hg38.bam"],
    "epiformat":"modbam",
    "cpg_coordinates": "tests/data/modbam_test-fw-rev.cpg_file.0based.hg38.bed.gz"
    }
    config_dict.update(test_specific)
    runner = EpiToBedgraph(config_dict)
    runner.read_mixture()
    runner.calc_coverage()
    #assert runner.coverage[0][0]==11
    assert (sum(runner.coverage[0]) == 4) and (sum(runner.coverage[1]) == 4)
    
