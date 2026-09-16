###################################################
#
# Script: ModbamReader.py
# Author: Ben Berman
# Description: parse epiread files for EM
#
# Copyright (c) 2023 ben berman, base4analytics, and volitionRx
# ###################################################

import pysam
from epiread_tools.epiparser import EpireadReader
from  epiread_tools.em_utils import cut
from epiread_tools.naming_conventions import methylation_state
import scipy.sparse as sp
import numpy as np


class ModbamReader(EpireadReader):

    class_epireads_processed = 0;
    class_mods_processed = 0;

    def __init__(self, fp):
        super().__init__(fp)
#        self.row = PatRow
        self.mapper_slop = 500 # Does this actually do anything?  CpGs outside the markers are not included in deconvolution anyway.
        self.min_qual_score = 170; # 230; # This is out of 255.  qual=256*prob. q(230)=0.9. q(170)=0.667 Probably need something more complicated if we use C/hmC/mC model

    def to_csr(self, epi_file, mapper):

        row = []
        col = []
        data = []

        epiread_iterator = cut(epi_file,
                               mapper.merge_intervals(mapper.original_intervals),
                               "bam")  # should be merged, unslopped
        row_ind = 0


        # pysam.AlignedSegment.modified_bases
        # Modified bases annotations from Ml/Mm tags. The output is 
        # Dict[(canonical base, strand, modification)] -> [ (pos,qual), …] 
        # with qual being (256*probability), or -1 if unknown. 
        # Strand==0 for forward and 1 for reverse strand modification



        for i, epiread in enumerate(epiread_iterator):
            ModbamReader.class_epireads_processed += 1

            mods = epiread.modified_bases # If you need coordinates relative to reverse_strand reads, *** use epiread.modified_bases_forward ***
#            print(f"epiread: {epiread.qname} mods: {mods}")


            if not mods:
                # TODO: *** I'm not sure if this is the right thing to do.  Check why the read has no mods
#                print(f"WARNING: epiread: {epiread.qname} has no mods")
                dummy = 1
            else:
                ref_positions = epiread.get_reference_positions(full_length=True) # *** I am not sure if full_length is right. It would make a huge difference for softclip, insertions.
#                print(f"reference_positions (len={len(ref_positions)}): {ref_positions}")
                # With full-length=True, we get some "None" values at the ends for softclipped bases
                #17732441, 17732442, 17732443, 17732444, 17732445, 17732446, 17732447, 17732448, 17732449, 17732450, 17732451, 17732452, 17732453, 17732454, 17732455, 17732456, 17732457, 17732458, 17732459, 17732460, 17732461, 17732462, 17732463, 17732464, 17732465, 17732466, 17732467, 17732468, 17732469, 17732470, 17732471, 17732472, 17732473, 17732474, 17732475, 17732476, 17732477, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None

                for key, values in mods.items():
#                    print(f"read: {epiread.qname} key: {key}\tvalues: {values}")
                    rev_strand = key[1] # 0 for forward, 1 for reverse
    #            for values in mods.values(): # I think we take the forward and reverse strands together
                    for pos, qual in values:
                        # Process pos and qual here
#                        print(f"Position: {pos}\tQual: {qual}")

                        # Filter quality
                        # TODO ***  I'm not sure if this logic is correct, in the C, 5mC, 5hmC model.  Check modkit
                        meth = methylation_state["N"]
                        if (qual > self.min_qual_score):
                            meth = methylation_state["C"]
                        elif (qual < (256 - self.min_qual_score)):
                            meth = methylation_state["T"]

                        # Get the absolute position of the base
                        abs = ref_positions[pos] # This is the reference position of the base
#                        print(f"Recording modification {abs},q={qual}")

                        # We should never call a modification on a base that is not in the reference?
                        # Actually, soft clips appear as None so we can skip these
                            # assert abs != None 
                        if (abs != None):
                            ModbamReader.class_mods_processed += 1

                            if rev_strand:
                                abs = abs - 1 # This is the way it is done for epireads

                            rel = mapper.abs_to_rel.get(abs)
                            if (rel):
                                if rel in mapper.all_rel:
                                    row.append(row_ind)  # Append the current row_ind
                                    col.append(mapper.abs_to_ind(abs))
                                    data.append(meth)
                                    row_ind+=1  # Increment row_ind after appending the rows

        # Last thing
        print (f"\tglobal epireads processed: {ModbamReader.class_epireads_processed}, global mods processed: {ModbamReader.class_mods_processed}")
                                    

        if not len(data): #no coverage
            return None, None

        return sp.csr_matrix((data, (row, col)), shape=(row_ind + 1, mapper.max_cpgs), dtype=int), np.zeros(shape=row_ind+1)

    
