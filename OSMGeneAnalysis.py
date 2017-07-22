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

import numpy as np

from collections import namedtuple
import sys

from OSMGeneEvidence import GenomeEvidence


class GeneDictionary(object):

    def __init__(self, args, log, genome_gff):
        # Shallow copies of the runtime environment.
        self.log = log
        self.args = args
        self.ContigRecord = namedtuple("ContigRecord", "Contig Genedict Genelist")
        self.contig_dict = self.genome_contig_dict(genome_gff)

    def get_gene_id(self, contig_id, gene_id):

        if contig_id in self.contig_dict:
            if gene_id in self.contig_dict[contig_id].Genedict:
                return self.contig_dict[contig_id].Genedict[gene_id]
        return None

    def genome_contig_dict(self, genome_gff):

        contig_dict = {}
        for contig_seqrecord in genome_gff:

            gene_list = self.sorted_gene_list(contig_seqrecord)  # sorted by increasing contig sequence position
            gene_dict = {}
            for gene in contig_seqrecord.features:  # process all the features in a biopython seqrecord contig object.
                gene_dict[gene.id] = gene

            contig_dict[contig_seqrecord.id] = self.ContigRecord(*[contig_seqrecord, gene_dict, gene_list])

        return contig_dict

    def sorted_gene_list(self, contig_seqrecord):

        gene_list = []
        for gene in contig_seqrecord.features:
            gene_list.append(gene)

        return sorted(gene_list, key=lambda g: g.location.start)

    def print_contig_dict(self):

        for contigkey, contigvalue in self.contig_dict.items():
            print("\n******* Contig:", contigvalue.Contig)
            for gene in contigvalue.Genelist:
                print("\n%%% Gene:", gene)
