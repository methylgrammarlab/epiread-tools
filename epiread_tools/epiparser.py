###################################################
#
# Script: epiparser.py
# Author: Irene Unterman
# Description: parse epiread files for EM
###################################################

import numpy as np
import pysam
from  epiread_tools.naming_conventions import *
from  epiread_tools.em_utils import Mapper, find_intersection, GenomicInterval
import scipy.sparse as sp
import os
import re
# import sys
# sys.path.insert(0, os.path.abspath('..'))
#%%

def tabix_verify(fp):
    print(os.getcwd())
    assert os.path.isfile(fp)  # file exists
    assert fp.endswith(".gz")  # gzipped file
    assert os.path.isfile(fp + ".tbi")  # index file exists

class Parser:
    '''
    read epiread data into scipy matrix
    each row is read even if it is empty e.g. NN-N
    '''

    def __init__(self,chrom, intervals, epiread_files, CpG_file, epi_format):
        '''

        :param chrom: chromosome name
        :param intervals: GenpmicInterval list
        :param epiread_files: list of filepaths
        :param CpG_file: file with coordinates of cpgs in genome
        :param epi_format: epiread file format
        '''
        self.chrom = chrom
        self.intervals = intervals
        self.epiread_files = epiread_files
        tabix_verify(CpG_file)
        self.CpG_file = CpG_file
        self.init_mapper()
        self.epiformat = epi_format
        self.fileobj= format_to_fileobj[self.epiformat]

    def init_mapper(self):
        '''
        map matrix coordinates
        :return:
        '''
        self.mapper = Mapper(self.chrom, self.intervals, self.epiread_files, self.CpG_file)

    def parse(self):
        '''
        read files at designated intervals to matrix
        :return: scipy sparse matrix, mapper
        '''
        small_matrices = []
        sources = []
        i = 0
        for epi_file in self.epiread_files:
            small_matrix = self.fileobj(epi_file).to_csr(self.mapper, self.intervals)
            small_matrices.append(small_matrix)
            sources.append((i, i + small_matrix.shape[0], self.mapper.sample_to_id[epi_file]))
            i += small_matrix.shape[0]
        self.mapper.init_index_to_source(sources)
        methylation_matrix = sp.vstack(small_matrices, dtype=int)
        return methylation_matrix, self.mapper
#%%
#epiread objects

class Epiread_format:
    '''
    epiread format includes: Chromosome name, Read name
    Read position in paired-end sequencing, Bisulfite strand,
    Position of the cytosine in the first CpG (0-based),  Retention pattern,
    Position of the first SNP, Base call of all SNPs covered
    min_start  is min(firstSNP,firstCpG) and max_end is max(lastSNP,lastCpG)
    '''
    def __init__(self, fp):
        tabix_verify(fp)
        self.fp = fp
        self.row = EpiRow

    def cut(self, intervals):
        '''
        Tabix out interval from file, return reads
        in original format
        :param intervals: list of GenomicInterval of chunk
        :return: iterator of relevant records
        '''
        tabixfile = pysam.TabixFile(self.fp)
        for interval in intervals:
            try:
                yield from tabixfile.fetch(interval.chrom, interval.start, interval.end)
            except ValueError:  # no coverage for this region in file
                continue

    def to_csr(self, mapper, intervals):
        row = []
        col = []
        data = []
        i = 0
        rel_intervals = np.array(mapper.rel_intervals).flatten()
        epiread_iterator = self.cut(intervals)
        for i, epiread in enumerate(epiread_iterator):
            record = self.row(*epiread.split(TAB))
            rel_start = mapper.abs_to_rel[record.get_start()]
            rel_end = rel_start+len(record)
            for intersect_start, intersect_end in find_intersection(rel_intervals, rel_start, rel_end):
                if intersect_end - intersect_start:  # overlap length is > 0
                    row.extend([i] * (intersect_end - intersect_start))
                    col.extend(list(range(mapper.rel_to_ind(intersect_start), mapper.rel_to_ind(intersect_end - 1) + 1)))
                    for cpg in record.methylation[intersect_start - rel_start:intersect_end - rel_start]:
                        data.append(methylation_state[cpg])
        return sp.csr_matrix((data, (row, col)), shape=(i + 1, mapper.max_cpgs), dtype=int)

class CoordsEpiread(Epiread_format):

    def __init__(self, fp):
        super().__init__(fp)
        self.row = CoordsRow

    def to_csr(self, mapper, intervals):
        row = []
        col = []
        data = []
        epiread_iterator = self.cut(intervals)
        i=0
        for i, epiread in enumerate(epiread_iterator):
            record = self.row(*epiread.split(TAB))
            for abs, cpg in record.get_coord_methylation():
                if mapper.abs_to_rel[abs] in mapper.all_rel:
                    row.append(i)
                    col.append(mapper.abs_to_ind(abs))
                    data.append(methylation_state[cpg])
        return sp.csr_matrix((data, (row, col)), shape=(i + 1, mapper.max_cpgs), dtype=int)

class EpiSNP(Epiread_format):

    def __init__(self, fp):
        super().__init__(fp)
        self.row = SNPRow

    def to_csr(self, mapper, intervals): #problem: aligning
        epiread_iterator = self.cut(intervals)
        row = []
        col = []
        data = []
        i = 0
        for i, epiread in enumerate(epiread_iterator):
            record = self.row(*epiread.split(TAB))
            if record.has_snps():
                rel_start = mapper.snp_abs_to_rel(record.get_snp_start())
                for dist, snp in record.get_snps():
                    if rel_start + int(dist) < mapper.max_snps and rel_start + int(dist):  # doesn't pass interval
                        row.append(i)
                        col.append(rel_start + int(dist))
                        data.append(dna_to_ind[snp])
        return sp.csr_matrix((data, (row, col)), shape=(i+1, mapper.max_snps), dtype=int)

class CommentEpiSNP(EpiSNP):
    def __init__(self, fp):
        super().__init__(fp)
        self.row = CommentEpiSNPRow

#%%
#Row objects

class EpiRow:
    '''
    epiread format includes: Chromosome name, Read name
    Read position in paired-end sequencing, Bisulfite strand,
    Position of the cytosine in the first CpG (0-based),  Retention pattern,
    Position of the first SNP, Base call of all SNPs covered
    min_start  is min(firstSNP,firstCpG) and max_end is max(lastSNP,lastCpG)
    '''
    def __init__(self, chrom, min_start, max_end, read_name, read_pos, strand, read_start, methylation,
                 snp_start="", snps="", origin=""):
        self.chrom = chrom
        self.min_start = min_start
        self.max_end = max_end
        self.read_name = read_name
        self.read_pos = read_pos
        self.strand = strand
        self.read_start = read_start
        self.methylation = methylation
        self.snp_start = snp_start
        self.snps = snps
        self.origin=origin

    def get_start(self):
        '''
        :return: read start e.g. 2100056
        '''
        return int(self.read_start)

    def has_snps(self):
        '''
        :return: True if SNPs in read
        '''
        #check that snp_start isn't "."
        return len(self.snp_start) and (self.snp_start != NO_DATA)

    def __len__(self):
        return len(self.methylation)

    def __repr__(self):
        return (TAB).join([self.chrom, self.read_name, str(self.read_pos), self.strand,\
        str(self.read_start), self.methylation, str(self.snp_start),\
        self.snps])

class CoordsRow(EpiRow):
    def __init__(self, chrom, min_start, max_end, read_name, read_pos, strand, coords, methylation,
                 snp_start="", snps="", origin=""):
        self.coords = [int(x) for x in coords.split(COORD_SEP)]
        self.read_start = self.coords[0]
        super().__init__(chrom, min_start, max_end, read_name, read_pos, strand, self.read_start, methylation,
                 snp_start, snps, origin)
        assert (len(self.coords) == len(self.methylation),
                "unequal length of coordinates and methylation, did you mean old_epiread?")

    def get_coord_methylation(self):
        '''
        pairs of coordinate, methylation
        e.g. 2051767,C
        :return:
        '''
        yield from zip(self.coords, self.methylation)

    def get_end(self):
        '''
        get last coordinate
        :return:
        '''
        return self.coords[-1]

class SNPRow(EpiRow):

    def __init__(self, chrom, min_start, max_end, read_name, read_pos, strand, read_start, methylation,
                 snp_start="", snps="", origin=""):
        super().__init__(chrom, min_start, max_end, read_name, read_pos, strand, read_start, methylation,
                 snp_start, snps, origin)
        self.snps = snps.split(SNP_SEP)

    def get_snps(self):
        '''
        :return: pairs of dist from start, SNP
        e.g. 0:A
        '''
        for i in range(0, len(self.snps), 2):
            yield (self.snps[i],self.snps[i+1])

    def get_snp_start(self):
        '''
        :return: snp start, e.g. 213356
        '''
        try:
            return int(self.snp_start)
        except ValueError:
            raise("Record has no SNPs")


class CommentEpiSNPRow(SNPRow):
    def __init__(self, chrom, min_start, max_end, read_name, read_pos, strand, read_start, methylation,
                 snp_start="", snps="", origin=""):
        super().__init__(chrom, min_start, max_end, read_name, read_pos, strand, read_start, methylation,
                 snp_start, snps, origin)
        self.snps = re.sub(r"[\(\[].*?[\)\]]", "", snps).split(SNP_SEP)

format_to_fileobj = {"old_epiread":Epiread_format, "old_epiread_A": CoordsEpiread, "clean_snps": EpiSNP,
                    "snps_with_comments": CommentEpiSNP
                     }

