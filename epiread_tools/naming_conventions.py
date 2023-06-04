#############################################################
# FILE: naming_conventions.py
# WRITER: Irene Unterman
# DESCRIPTION: Contains all conventions for using names
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

#%%
from collections import defaultdict
initial_high, initial_low = 0.8, 0.15 #initialization for EM
pseudocount = 1e-10 #same as methpipe
NOVAL, UNMETHYLATED, METHYLATED = range(3)
A,C,T,G, R, Y = range(1, 7)
dna_to_ind = defaultdict(lambda: 0)
dna_to_ind["A"], dna_to_ind["C"], dna_to_ind["T"], dna_to_ind["G"], \
dna_to_ind["R"],dna_to_ind["Y"] = A, C, T, G, R, Y
#R = C/T on positive strand
#Y = G/A on negative strand
ind_to_dna = {NOVAL: ".", A:"A", C:"C", T:"T", G:"G", R: "R", Y: "Y"}
methylation_state = {"C": METHYLATED, "T": UNMETHYLATED, "c": METHYLATED,
                     "t": UNMETHYLATED, ".": NOVAL, "N": NOVAL, "-":NOVAL, "n": NOVAL}
NO_DATA = "."
SNP_SEP = ":"
COORD_SEP = ","
POSITIVE_STRAND = "+"
NEGATIVE_STRAND = "-"
TAB = "\t"
ALL = 0
upper_conf_thresh = 0.9
lower_conf_thresh = 0.1


CHROMOSOMES = ["chr" + str(x) for x in [*range(1,23), "X", "Y", "M"]]
atlas_formats = ["meth_cov", "beta_matrices", "lambda_matrices"]