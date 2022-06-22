#############################################################
# FILE: naming_conventions.py
# WRITER: Irene Unterman
# DESCRIPTION: Contains all conventions for using names
#############################################################
#%%
from collections import defaultdict
initial_high, initial_low = 0.8, 0.15
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