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

from __future__ import print_function,  division

from collections import namedtuple

from OSMGeneEvidence import GenomeEvidence


class GeneDictionary(object):

    def __init__(self, args, log, genome_gff):
        # Shallow copies of the runtime environment.
        self.log = log
        self.args = args
        self.ContigRecord = namedtuple("ContigRecord", "Contig Genedict Genelist CDSlist")
        self.contig_dict = self.genome_contig_dict(genome_gff)

    def genome_contig_dict(self, genome_gff):

        contig_dict = {}
        for contig_seqrecord in genome_gff:

            gene_dict = self.get_gene_dict(contig_seqrecord)  # lookup using gene id.
            gene_list = self.sorted_gene_list(contig_seqrecord)  # sorted by increasing contig sequence position
            CDS_list = self.sorted_CDS_list(contig_seqrecord)  # sorted by increasing contig sequence position

            contig_dict[contig_seqrecord.id] = self.ContigRecord(*[contig_seqrecord, gene_dict, gene_list, CDS_list])

        return contig_dict

    def get_gene_dict(self,  contig_seqrecord):

        gene_dict = {}
        for gene in contig_seqrecord.features:  # process all the features in a biopython seqrecord contig object.
            gene_dict[gene.id] = gene

        return gene_dict

    def sorted_gene_list(self, contig_seqrecord):

        gene_list = []
        for gene in contig_seqrecord.features:
            gene_list.append(gene)

        return sorted(gene_list, key=lambda g: g.location.start)

    def sorted_CDS_list(self, contig_seqrecord):

        CDS_list = []
        for gene in contig_seqrecord.features:
            CDS_list = CDS_list + self.get_CDS_list(gene.sub_features)

        return sorted(CDS_list, key=lambda cds: cds.location.start)

    def get_CDS_list(self, sub_feature_list):  # recursively descend the feature tree and create a list of CDS features.

        CDS_list = []

        if sub_feature_list is not None:

            for sub_feature in sub_feature_list:

                if sub_feature.type == "CDS":

                    CDS_list.append(sub_feature)

                CDS_list = CDS_list + self.get_CDS_list(sub_feature.sub_features)

        return CDS_list

    def print_contig_dict(self):

        for contigkey, contigvalue in self.contig_dict.items():
            print("\n******* Contig:", contigvalue.Contig)
            for gene in contigvalue.CDSlist:
                print("\n%%% CDS:", gene)


class GeneAnalysis(object):

    def __init__(self, args, log, genome_gff, sam_file_name):
        # Shallow copies of the runtime environment.
        self.log = log
        self.args = args
        self.gene_dictionary = GeneDictionary(self.args, self.log, genome_gff)
        self.gene_evidence = GenomeEvidence(self.args, self.log, genome_gff)
        self.gene_evidence.read_sam_file(sam_file_name)

