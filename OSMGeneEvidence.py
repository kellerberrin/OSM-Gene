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
import copy

from OSM_Filter import SNPAnalysis

######################################################################################################################
# Object to return the SAM file read evidence.
# Nucleotide counts and deletions go into an unsigned integer numpy area the same size as the contig. region sequence.
# Insertions are entered into a numpy array of objects initially set to None. Insertions are added at the
# appropriate sequence index as dictionaries indexed by inserted sequence with count as value.
######################################################################################################################


class GenomeEvidence(object):

    nucleotide_offset = {"A": 0, "a" : 0, "C": 1, "c": 1, "G": 2, "g" : 2, "T": 3, "t": 3, "U": 3, "u" : 3, "-": 4, "+": 5}
    nucleotide_list = [ "A", "C", "G", "T", "-", "+"]
    EvidenceFields = namedtuple("EvidenceFields", "Contigrecord Contigfixedarray Contiginsertarray")

    def __init__(self, log, sam_evidence):

        self.log = log
        self.genome_evidence = self.__create_genome_evidence(sam_evidence)


    ##############################################################################################
    #
    #   Public class members
    #
    ##############################################################################################

    def get_all_snp(self, min_read_count, min_mutant_proportion):
        snp_evidence = self.__generate_all_snp(min_read_count, min_mutant_proportion)
        return SNPAnalysis(self.log, self.genome_evidence, snp_evidence)

    ##############################################################################################
    #
    #   Private class members
    #
    ##############################################################################################

    def __generate_all_snp(self, min_read_count, min_mutant_proportion):

        self.log.info("Filtering: %d contiguous regions for SNP, minimum reads: %d, minimum mutant proportion: %f"
                      , len(self.genome_evidence), min_read_count, min_mutant_proportion)

        snp_evidence = {}
        for contig_id, contig_evidence in self.genome_evidence.items():

            contig_snp_list = []
            fixed = contig_evidence.Contigfixedarray
            contig = contig_evidence.Contigrecord

            for idx in range(len(contig.seq)):

                nucleotide = contig.seq[idx]
                nucleotide_count = fixed[idx][GenomeEvidence.nucleotide_offset[nucleotide]]
                a_count = fixed[idx][GenomeEvidence.nucleotide_offset["A"]]
                c_count = fixed[idx][GenomeEvidence.nucleotide_offset["C"]]
                g_count = fixed[idx][GenomeEvidence.nucleotide_offset["G"]]
                t_count = fixed[idx][GenomeEvidence.nucleotide_offset["T"]]  # nucleotide 'T' and 'U'
                delete_count = fixed[idx][GenomeEvidence.nucleotide_offset["-"]]
                insert_count = fixed[idx][GenomeEvidence.nucleotide_offset["+"]]
                count_list = [a_count, c_count, g_count, t_count, delete_count, insert_count]
                sum_count_list = float(sum(count_list))
                if sum_count_list > 0:
                    proportion_mutant = 1.0 - (nucleotide_count / sum_count_list)
                else:
                    proportion_mutant = 0.0

                if proportion_mutant >= min_mutant_proportion and sum_count_list >= min_read_count:
                    contig_snp_list.append(idx)

            snp_evidence[contig.id] = SNPAnalysis.SNPFields(Contigrecord=contig, SNPlist=contig_snp_list)
            self.log.info("Contig: %s has %d raw SNP locations", contig.id, len(contig_snp_list))

        return snp_evidence

    def __create_genome_evidence(self, sam_evidence):

        evidence = {}
        for contig_id, contig_evidence in sam_evidence.items():
            # The fixed array in sam_evidence is a c_type RawArray so we will do a deepcopy
            contig_fixed_array = copy.deepcopy(contig_evidence.Contigfixedarray)
            evidence[contig_id] = GenomeEvidence.EvidenceFields( Contigrecord=contig_evidence.Contigrecord
                                                                , Contigfixedarray=contig_fixed_array
                                                                , Contiginsertarray=contig_evidence.Contiginsertarray)
        return evidence

