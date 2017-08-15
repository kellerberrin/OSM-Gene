# MIT License
#
# Copyright (c) 2017 
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
#
from __future__ import absolute_import, division, print_function, unicode_literals

from collections import namedtuple
import bisect
import sys

from OSMGeneEvidence import GenomeEvidence


##############################################################################################
#
#   A contig indexed dictionary of Genes and CDS locations sorted by increasing sequence.
#
##############################################################################################


class GeneDictionary(object):

    def __init__(self, genome_evidence):

        self.GeneDict = namedtuple("GeneDict", "Contigrecord Geneoffsets Genelist Cdsoffsets Cdslist")
        self.contig_dict = self.__genome_contig_dict(genome_evidence.get_genome_evidence())

    ##############################################################################################
    #
    #   Public class members
    #
    ##############################################################################################

    def get_gene_dictionary(self):
        return self.contig_dict

    def get_gene(self, contig_id, sequence_idx):  # returns NULL if the sequence_idx does not lie within a gene.
        return self.__lookup_gene_feature_list(contig_id, sequence_idx)

    def get_cds(self, contig_id, sequence_idx):  # returns NULL if the sequence_idx does not lie within a cds
        return self.__lookup_cds_feature_list(contig_id, sequence_idx)

    ##############################################################################################
    #
    #   Private class members
    #
    ##############################################################################################

    def __genome_contig_dict(self, genome_evidence):

        contig_dict = {}
        for contig_id, contig_evidence in genome_evidence.items():
            contig_record = contig_evidence.Contigrecord
            gene_list = self.__sorted_gene_list(contig_record)  # sorted by increasing contig sequence position
            gene_offset_list = self.__get_feature_offsets(gene_list)  # sorted list of gene sequence index offsets.
            cds_list = self.__sorted_cds_list(contig_record)  # sorted by increasing contig sequence position
            cds_offset_list = self.__get_feature_offsets(cds_list)  # sorted list of gene sequence index offsets.

            contig_dict[contig_id] = self.GeneDict(Contigrecord=contig_record
                                                    , Geneoffsets=gene_offset_list
                                                    , Genelist=gene_list
                                                    , Cdsoffsets=cds_offset_list
                                                    , Cdslist=cds_list)

        return contig_dict


    def __sorted_gene_list(self, contig_seqrecord):

        gene_list = []
        for gene in contig_seqrecord.features:
            gene_list.append(gene)

        return sorted(gene_list, key=lambda g: g.location.start)

    def __get_feature_offsets(self, sorted_feature_list):

        feature_offset_list = []
        for feature in sorted_feature_list:  # process all the features in a biopython seqrecord contig object.
            feature_offset_list.append(feature.location.start)

        return feature_offset_list

    def __lookup_gene_feature_list(self, contig_id, sequence_idx):

        def find_le(a, x):
            'Find rightmost value less than or equal to x'
            i = bisect.bisect_right(a, x)
            return (i-1) if i else None

        gene_dict = self.contig_dict[contig_id]
        gene_idx = find_le(gene_dict.Geneoffsets, sequence_idx)
        if gene_idx is not None:
            gene = gene_dict.Genelist[gene_idx]
            if sequence_idx < gene.location.start or sequence_idx > gene.location.end:
                gene = None
        else:
            gene = None

        return gene

    def __lookup_cds_feature_list(self, contig_id, sequence_idx):

        def find_le(a, x):
            'Find rightmost value less than or equal to x'
            i = bisect.bisect_right(a, x)
            return (i-1) if i else None

        gene_dict = self.contig_dict[contig_id]
        cds_idx = find_le(gene_dict.Cdsoffsets, sequence_idx)
        if cds_idx is not None:
            cds = gene_dict.Cdslist[cds_idx]
            if sequence_idx < cds.location.start or sequence_idx > cds.location.end:
                cds = None
        else:
            cds = None

        return cds

    def __sorted_cds_list(self, contig_seqrecord):

        cds_list = []
        for gene in contig_seqrecord.features:
            cds_list = cds_list + GeneDictionary.__get_cds_list(gene.sub_features)

        return sorted(cds_list, key=lambda cds: cds.location.start)

    @staticmethod
    def __get_cds_list(sub_feature_list):  # recursively descend the feature tree and create a list of cds features.

        cds_list = []

        if sub_feature_list is not None:

            for sub_feature in sub_feature_list:

                if sub_feature.type == "cds" or sub_feature.type == "CDS":

                    cds_list.append(sub_feature)

                cds_list = cds_list + GeneDictionary.__get_cds_list(sub_feature.sub_features)

        return cds_list


##############################################################################################
#
#   Generates SNP evidence from the supplied genome_evidence.
#
##############################################################################################


class GenomeSNPFilter(object):

    def __init__(self, log, genome_evidence, min_read_count=20, min_mutant_proportion=0.7):

        self.log = log
        self.SNPFields = namedtuple("SNPFields", "Contigrecord SNPlist")
        self.genome_evidence = genome_evidence
        self.gene_dictionary = GeneDictionary(genome_evidence)
        self.snp_evidence = self.filter_count_proportion(min_read_count, min_mutant_proportion)

    def get_snp_evidence(self):
        return self.snp_evidence

    def set_snp_evidence(self, snp_evidence):
        self.snp_evidence = snp_evidence

    def filter_count_proportion(self, min_read_count, min_mutant_proportion):

        self.log.info("Filtering: %d contiguous regions for SNP, minimum reads: %d, minimum mutant proportion: %f"
                      , len(self.genome_evidence.get_genome_evidence()), min_read_count, min_mutant_proportion)

        snp_evidence = {}
        for contig_id, contig_evidence in self.genome_evidence.get_genome_evidence().items():

            contig_snp_list = []
            fixed = contig_evidence.Contigfixedarray
            insert = contig_evidence.Contiginsertarray
            contig = contig_evidence.Contigrecord

            for idx in range(len(contig.seq)):

                nucleotide = contig.seq[idx]
                nucleotide_count = fixed[idx][GenomeEvidence.nucleotide_offset[nucleotide]]
                a_count = fixed[idx][GenomeEvidence.nucleotide_offset["A"]]
                c_count = fixed[idx][GenomeEvidence.nucleotide_offset["C"]]
                g_count = fixed[idx][GenomeEvidence.nucleotide_offset["G"]]
                t_count = fixed[idx][GenomeEvidence.nucleotide_offset["T"]]  # nucleotide 'T' and 'U'
                delete_count = fixed[idx][GenomeEvidence.nucleotide_offset["-"]]
                insert_count = 0 if insert[idx] is None else len(insert[idx])
                count_list = [a_count, c_count, g_count, t_count, delete_count, insert_count]
                sum_count_list = float(sum(count_list))
                if sum_count_list > 0:
                    proportion_mutant = 1.0 - (nucleotide_count / sum_count_list)
                else:
                    proportion_mutant = 0.0

                if proportion_mutant >= min_mutant_proportion and sum_count_list >= min_read_count:
                    contig_snp_list.append(idx)

            snp_evidence[contig.id] = self.SNPFields(Contigrecord=contig, SNPlist=contig_snp_list)
            self.log.info("Contig: %s has %d raw SNP locations", contig.id, len(contig_snp_list))

        return snp_evidence

    def filter_gene_snp(self):

        self.log.info("Filtering: %d contiguous regions for SNP within gene boundaries"
                      , len(self.genome_evidence.get_genome_evidence()))

        snp_evidence = {}
        for contig_id, contig_snp in self.snp_evidence.items():

            gene_snp_list = []
            for snp in contig_snp.SNPlist:
                gene =self.gene_dictionary.get_gene(contig_id, snp)
                if gene is not None:
                    gene_snp_list.append(snp)
            self.log.info("Contig: %s has %d SNP within gene boundaries", contig_id, len(gene_snp_list))
            snp_evidence[contig_id] = self.SNPFields(Contigrecord=contig_snp.Contigrecord, SNPlist=gene_snp_list)

        return snp_evidence

    def filter_cds_snp(self):

        self.log.info("Filtering: %d contiguous regions for SNP within cds boundaries"
                      , len(self.genome_evidence.get_genome_evidence()))

        snp_evidence = {}
        for contig_id, contig_snp in self.snp_evidence.items():

            cds_snp_list = []
            for snp in contig_snp.SNPlist:
                cds =self.gene_dictionary.get_cds(contig_id, snp)
                if cds is not None:
                    cds_snp_list.append(snp)
            self.log.info("Contig: %s has %d SNP within cds boundaries", contig_id, len(cds_snp_list))
            snp_evidence[contig_id] = self.SNPFields(Contigrecord=contig_snp.Contigrecord, SNPlist=cds_snp_list)

        return snp_evidence