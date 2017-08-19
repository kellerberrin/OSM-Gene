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

import sys
from collections import namedtuple
import bisect


##############################################################################################
#
#   An auxiliary helper class to provide indexed look up on contiguous features.
#   A contiguous indexed dictionary of Genes and CDS locations sorted by increasing sequence.
#
##############################################################################################


class GeneDictionary(object):

    GeneDict = namedtuple("GeneDict", "Contigrecord Geneoffsets Genelist Cdsoffsets Cdslist")

    def __init__(self, genome_evidence):

        self.contig_dict = self.__genome_contig_dict(genome_evidence)

    ##############################################################################################
    #
    #   Public class members
    #
    ##############################################################################################

    def get_gene(self, contig_id, sequence_idx):  # returns NULL if the sequence_idx does not lie within a gene.
        return self.__lookup_gene_feature_list(contig_id, sequence_idx)

    def get_cds(self, contig_id, sequence_idx):  # returns NULL if the sequence_idx does not lie within a cds
        return self.__lookup_cds_feature_list(contig_id, sequence_idx)

    def get_upstream_gene(self, contig_id, sequence_idx):  # returns the nearest gene in increasing contig. sequence.
        return self.__upstream_feature_list(contig_id, sequence_idx)

    def get_downstream_gene(self, contig_id, sequence_idx):  # returns the nearest gene in decreasing contig sequence.
        return self.__downstream_feature_list(contig_id, sequence_idx)

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

            contig_dict[contig_id] = GeneDictionary.GeneDict(Contigrecord=contig_record
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

        return gene  # returns None if the sequence_idx is not within a gene, else returns the CDS (SeqFeature) record.

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

        return cds  # returns None if the sequence_idx is not within a CDS, else returns the CDS (SeqFeature) record.

    def __upstream_feature_list(self, contig_id, sequence_idx):

        def find_ge(a, x):
            'Find leftmost item greater than or equal to x'
            i = bisect.bisect_left(a, x)
            return i if i != len(a) else None

        gene_dict = self.contig_dict[contig_id]
        gene_idx = find_ge(gene_dict.Geneoffsets, sequence_idx)
        if gene_idx is not None:
            gene = gene_dict.Genelist[gene_idx]
        else:
            gene = None

        return gene  # returns None if the sequence_idx is more than the beginning of the last gene in the contig.


    def __downstream_feature_list(self, contig_id, sequence_idx):

        def find_le(a, x):
            'Find rightmost value less than or equal to x'
            i = bisect.bisect_right(a, x)
            return (i-1) if i else None

        gene_dict = self.contig_dict[contig_id]
        gene_idx = find_le(gene_dict.Geneoffsets, sequence_idx)
        if gene_idx is not None:
            gene = gene_dict.Genelist[gene_idx]
        else:
            gene = None

        return gene  # returns None if the sequence_idx is less than the beginning of the first gene in the contig.

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
#   Manipulates SNP evidence.
#
##############################################################################################


class SNPAnalysis(object):

    SNPFields = namedtuple("SNPFields", "Contigrecord SNPlist")

    def __init__(self, log, genome_evidence, snp_evidence):

        self.log = log
        self.genome_evidence = genome_evidence
        self.snp_evidence = snp_evidence
        self.gene_dictionary = GeneDictionary(genome_evidence)

    ##############################################################################################
    #
    #   Public class members
    #
    ##############################################################################################

    def get_snp_evidence(self):
        return self.snp_evidence

    def get_gene_dictionary(self):
        return self.gene_dictionary

    def get_genome_evidence(self):
        return self.genome_evidence

    def filter_gene_snp(self):

        self.log.info("Filtering: %d contiguous regions for SNP within gene boundaries", len(self.genome_evidence))
        snp_evidence = {}
        for contig_id, contig_snp in self.snp_evidence.items():

            gene_snp_list = []
            for snp in contig_snp.SNPlist:
                gene = self.gene_dictionary.get_gene(contig_id, snp)
                if gene is not None:
                    gene_snp_list.append(snp)
            self.log.info("Contig: %s has %d SNP within gene boundaries", contig_id, len(gene_snp_list))
            snp_evidence[contig_id] = SNPAnalysis.SNPFields(Contigrecord=contig_snp.Contigrecord
                                                            , SNPlist=gene_snp_list)

        return SNPAnalysis(self.log, self.genome_evidence,snp_evidence)

    def filter_cds_snp(self):

        self.log.info("Filtering: %d contiguous regions for SNP within cds boundaries", len(self.genome_evidence))
        snp_evidence = {}
        for contig_id, contig_snp in self.snp_evidence.items():

            cds_snp_list = []
            for snp in contig_snp.SNPlist:
                cds = self.gene_dictionary.get_cds(contig_id, snp)
                if cds is not None:
                    cds_snp_list.append(snp)
            self.log.info("Contig: %s has %d SNP within cds boundaries", contig_id, len(cds_snp_list))
            snp_evidence[contig_id] = SNPAnalysis.SNPFields(Contigrecord=contig_snp.Contigrecord
                                                            , SNPlist=cds_snp_list)

        return SNPAnalysis(self.log, self.genome_evidence,snp_evidence)

    def union(self, snp_analysis):
        def operator(s,t): return s.union(t)
        return self.__binary_set_operation(operator, snp_analysis)

    def intersection(self, snp_analysis):
        def operator(s,t): return s.intersection(t)
        return self.__binary_set_operation(operator, snp_analysis)

    def difference(self, snp_analysis):
        def operator(s,t): return s.difference(t)
        return self.__binary_set_operation(operator, snp_analysis)

    def symmetric_difference(self, snp_analysis):
        def operator(s,t): return s.symmetric_difference(t)
        return self.__binary_set_operation(operator, snp_analysis)

    ##############################################################################################
    #
    #   Private class members
    #
    ##############################################################################################

    def __binary_set_operation(self, operator, snp_analysis):
        snp_evidence = {}
        for contig_id, contig_snp in self.snp_evidence.items():

            if contig_id not in snp_analysis.snp_evidence:
                self.log.error("SNP union; contig region: %s not found in both snp evidence objects", contig_id)
                sys.exit()

            s_snp = set(contig_snp.SNPlist)
            t_snp = set(snp_analysis.snp_evidence[contig_id].SNPlist)
            result_snp_set = operator(s_snp, t_snp)
            result_snp_list = list(result_snp_set)
            snp_evidence[contig_id] = SNPAnalysis.SNPFields(Contigrecord=contig_snp.Contigrecord
                                                            , SNPlist=result_snp_list)
        return SNPAnalysis(self.log, self.genome_evidence,snp_evidence)
