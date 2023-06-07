###################################################
#
# Script: modbam_parser.py
# Author: Ben Berman
# Description: parse epiread files for EM
#
# Copyright (c) 2023 ben berman, base4analytics, and volitionRx
# ###################################################

import pysam
from epiread_tools.epiparser import EpiRow

# pysam.AlignedSegment.modified_bases
# Modified bases annotations from Ml/Mm tags. The output is 
# Dict[(canonical base, strand, modification)] -> [ (pos,qual), …] 
# with qual being (256*probability), or -1 if unknown. 
# Strand==0 for forward and 1 for reverse strand modification
