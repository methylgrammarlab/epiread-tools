###################################################
#
# Script: epiparser.py
# Author: Irene Unterman
# Description: parse epiread files for EM
###################################################

import numpy as np
import pandas as pd
import pysam
from  epiread_tools.naming_conventions import *
from  epiread_tools.em_utils import Mapper, find_intersection, bedgraph_to_intervals, cut
from  epiread_tools.em_utils import GenomicInterval, split_intervals_to_chromosomes, in_intervals
import scipy.sparse as sp
import os
import re

#%%

class EpireadReader:
    '''
    read epiread data into scipy matrix
    each row is read even if it is empty e.g. NN-N
    '''
    def __init__(self, config):
        self.config = config
        if config['bedfile']:
            self.genomic_intervals = bedgraph_to_intervals(config['genomic_intervals'], config['header'])
        else:
            self.genomic_intervals = [GenomicInterval(x) for x in config['genomic_intervals']]
        self.intervals_per_chrom = split_intervals_to_chromosomes(self.genomic_intervals)
        self.interval_order, self.matrices, self.cpgs = [],[],[]
        self.row = EpiRow

    def get_matrices_for_intervals(self):
        self.parse_multiple_chromosomes()
        return self.interval_order, self.matrices, self.cpgs

    def parse_multiple_chromosomes(self):
        for chrom, intervals in self.intervals_per_chrom.items():
            self.parse_one_chromosome(chrom, intervals)

    def parse_one_chromosome(self, chrom, intervals):
        '''
        get data amd mapper
        :return:
        '''
        intervals = sorted(intervals, key=lambda x: x.start)
        self.interval_order.extend(intervals)
        methylation_matrix, mapper = self.file_list_to_csr(chrom, intervals)
        window_list = mapper.get_ind_intervals(intervals)
        for start, end in window_list:
            slice = methylation_matrix[:,start:end]
            self.matrices.append(slice[slice.getnnz(1)>0]) #remove empty rows
            self.cpgs.append(np.array([mapper.ind_to_abs(x) for x in range(start, end)]))

    def file_list_to_csr(self, chrom, intervals):
        '''
        read files at designated intervals to matrix
        :return: scipy sparse matrix, mapper
        '''
        mapper = Mapper(chrom, intervals, self.config['epiread_files'], self.config['cpg_coordinates'])
        small_matrices = []
        sources = []
        i = 0
        for epi_file in self.config['epiread_files']:
            small_matrix = self.to_csr(epi_file, mapper)
            small_matrices.append(small_matrix)
            sources.append((i, i + small_matrix.shape[0], mapper.sample_to_id[epi_file]))
            i += small_matrix.shape[0]
        mapper.init_index_to_source(sources)
        methylation_matrix = sp.vstack(small_matrices, dtype=int)
        return methylation_matrix, mapper

    def to_csr(self, epi_file, mapper):
        row = []
        col = []
        data = []
        i = 0
        rel_intervals = np.array(mapper.rel_intervals).flatten()
        epiread_iterator = cut(epi_file, mapper.merged_intervals)
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


atlas_formats = ["meth_cov", "beta_matrices", "lambda_matrices"] #TODO: delete or move
'''
a note on strandedness:
All positions must match the cpg file. To process 
strand separately, make sure the CoG file includes separate coordinates
e.g. chrN:200-201, chrN:201-202. This has to be the same file the epireads 
are aligned to
'''
class AtlasReader:
    def __init__(self, config):
        self.config = config
        self.atlas = self.config['atlas_file']
        self.genomic_intervals = bedgraph_to_intervals(config['genomic_intervals'], config['header'])
        self.intervals_per_chrom = split_intervals_to_chromosomes(self.genomic_intervals)

    def meth_cov_to_beta_matrices(self):
        '''
        input for celfie+, list of matrices
        with beta values per region in regions
        :return: lsit of matrices
        '''
        self.load_meth_cov_atlas()
        interval_order, matrices = self.parse_multiple_chromosomes(self.beta)
        return matrices

    def meth_cov_to_sum(self):
        '''
        input for celfie, sum meth
        and cov for each TIM
        :return: sum meth per tim, sum cov per tim
        '''
        self.load_meth_cov_atlas()
        meth_interval_order, meth = self.parse_multiple_chromosomes(self.meth)
        cov_interval_order, cov = self.parse_multiple_chromosomes(self.cov)
        assert (meth_interval_order == cov_interval_order)
        sum_meth, sum_cov = np.array([np.nansum(x, axis=1) for x in meth]), np.array([np.nansum(x, axis=1) for x in cov])
        return sum_meth, sum_cov

    def load_meth_cov_atlas(self):
        '''
        read file with 3 BED cols, and meth+cov per cell type
        :return:
        '''
        df = pd.read_csv(self.atlas, sep="\t")
        vals = df.iloc[:, 3:].values
        self.meth = vals[:, ::2].T
        self.cov =  vals[:, 1::2].T
        self.beta = self.meth/self.cov
        self.atlas_chrom = df["CHROM"].values
        self.atlas_start = df["START"].values  # remove if 0 based

    def align_vals(self, chrom, mapper, vals):
        '''
        make sure cpgs are in intervals
        :param chrom: chromosome e.g. chr2
        :param mapper: Mapper instance
        :param vals: meth/cov/beta matrix
        :return: aligned matrix
        '''
        mat = np.ones(shape=(self.beta.shape[0], mapper.max_cpgs))  # cell types by cpgs
        mat.fill(np.nan)
        chrom_filter = (self.atlas_chrom == chrom)
        cpgs = self.atlas_start[chrom_filter]
        for i, cpg in enumerate(cpgs):
            if cpg in mapper.abs_to_rel and mapper.abs_to_rel[cpg] in mapper.all_rel:
                mat[:, mapper.abs_to_ind(cpg)] = vals[:, chrom_filter][:, i]
            elif not in_intervals(cpg, mapper.merged_intervals):
                print("not in intervals", cpg)
                pass
            else:
                raise KeyError(cpg, "atlas start doesn't match cpg file")
        return mat

    def parse_one_chromosome(self, chrom, intervals, vals):
        '''
        save vals per requested interval
        :param intervals: GenomicInterval ,may overlap
        :param chrom: chromosome e.g. chr2
        :param vals: meth/cov/beta matrix
        :return: list of matrix per interval
        '''
        matrices = []
        interval_order = []
        intervals = sorted(intervals, key=lambda x: x.start)
        interval_order.extend([str(x) for x in intervals])
        mapper = Mapper(chrom, intervals, [], self.config['cpg_coordinates'], False)  # init mapping
        window_list = mapper.get_ind_intervals(intervals)
        mat = self.align_vals(chrom, mapper, vals)
        for start, end in window_list:
            slice = mat[:, start:end]
            matrices.append(slice)
        return interval_order, matrices

    def parse_multiple_chromosomes(self, val):
        '''
        iterate over all chromosomes in intervals
        :param val:  meth/cov/beta matrix
        :return: interval order, list of matrix per interval
        '''
        interval_order = []
        matrices = []
        for chrom, intervals in self.intervals_per_chrom.items():
            sorted_intervals, mats = self.parse_one_chromosome(chrom, intervals, val)
            interval_order.extend(sorted_intervals)
            matrices.extend(mats)
        return interval_order, matrices


#%%
#epiread objects

class CoordsEpiread(EpireadReader):

    def __init__(self, fp):
        super().__init__(fp)
        self.row = CoordsRow

    def to_csr(self, epi_file, mapper):
        row = []
        col = []
        data = []
        epiread_iterator = cut(epi_file, mapper.merged_intervals)
        i=0
        for i, epiread in enumerate(epiread_iterator):
            record = self.row(*epiread.split(TAB))
            for abs, cpg in record.get_coord_methylation():
                if mapper.abs_to_rel[abs] in mapper.all_rel:
                    row.append(i)
                    col.append(mapper.abs_to_ind(abs))
                    data.append(methylation_state[cpg])
        return sp.csr_matrix((data, (row, col)), shape=(i + 1, mapper.max_cpgs), dtype=int)

class EpiSNP(EpireadReader):

    def __init__(self, fp):
        super().__init__(fp)
        self.row = SNPRow

    def to_csr(self,epi_file, mapper): #problem: aligning
        epiread_iterator = cut(epi_file, mapper.merged_intervals)
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


def tabix_verify(fp):
    assert os.path.isfile(fp)  # file exists
    assert fp.endswith(".gz")  # gzipped file
    assert os.path.isfile(fp + ".tbi")  # index file exists

epiformat_to_reader = {"old_epiread": EpireadReader, "old_epiread_A": CoordsEpiread}

#%%
# config = {"epiread_files": ["/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/EM_regions_100_6_rep9_mixture.epiread.gz"],
#           "atlas_file": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/EM_regions_100_atlas_over_tims.txt",
#           "genomic_intervals": "/Users/ireneu/PycharmProjects/deconvolution_models/tests/data/EM_regions_100_processed_tims.txt",
#           "cpg_coordinates":  "/Users/ireneu/PycharmProjects/in-silico_deconvolution/debugging/hg19.CpG.bed.sorted.gz",
#            "outfile":"sample_output.something",
#             "epiformat" : "old_epiread",
#           "bedfile":True,
#           "header":True,
#           "max_iterations":1000,
#           "random_restarts":1,
#           "stop_criterion":0.001
#           }
# em_model = AtlasReader(config)
# em_model.meth_cov_to_beta_matrices()