# ###################################################
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

import numpy as np
import pytest
from epiread_tools.naming_conventions import *
from epiread_tools.em_utils import GenomicInterval
from epiread_tools.epireadToBedgraph import EpiToBedgraph
from epiread_tools.epiparser import EpireadReader, CoordsEpiread, CommentEpiSNP

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

def test_bedgraph_from_intervals(config_dict):
    test_specific = {"genomic_intervals":["chr1:205499880-205500150","chr10:17599401-17600087"],
    "epiread_files" : ["tests/data/old_epiread.epiread.gz"],
                     "epiformat":"old_epiread"}
    config_dict.update(test_specific)
    runner = EpiToBedgraph(config_dict)
    runner.read_mixture()
    assert runner.matrices[0].shape == (10,3)

def test_A_option(config_dict):
    test_specific = {"genomic_intervals":["chr1:205499880-205500150","chr10:17599401-17600087"],
    "epiread_files" : ["tests/data/old_epiread_A.epiread.gz"],
                     "epiformat":"old_epiread_A"}
    config_dict.update(test_specific)
    runner = EpiToBedgraph(config_dict)
    runner.read_mixture()
    assert runner.matrices[0].shape == (8,3)


def test_overlapping_regions(config_dict):
    test_specific = {"genomic_intervals" : ["chr1:205499880-205500150", "chr1:205499880-205500150"],
    "epiread_files" : ["tests/data/old_epiread_A.epiread.gz"],
                     "epiformat":"old_epiread_A"}
    config_dict.update(test_specific)
    runner = EpiToBedgraph(config_dict)
    runner.read_mixture()
    runner.calc_coverage()
    assert all(runner.coverage[0]==runner.coverage[1])
    assert runner.coverage[0]==7




