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
from OSM_Filter import GeneDictionary, SNPAnalysis
from OSMGeneEvidence import GenomeEvidence




##############################################################################################
#
#   Compares the genome evidence to the fasta sequence or other genome evidence objects
#
##############################################################################################


class GeneAnalysis(object):

    def __init__(self, log):
        # Shallow copies of the runtime environment.
        self.log = log

    def print_snp_stats(self, snp_analysis_obj):

        total_SNP_mutations = 0
        for contig_id, contig_evidence in snp_analysis_obj.get_snp_evidence().items():

            self.log.info("********** Contig: %s, SNP: %d, *********", contig_id, len(contig_evidence.SNPlist))
            for snp in contig_evidence.SNPlist:
                self.print_snp(contig_id, snp, snp_analysis_obj.get_genome_evidence())
                gene = snp_analysis_obj.get_gene_dictionary().get_gene(contig_id, snp)
                if gene is not None:
                    cds = snp_analysis_obj.get_gene_dictionary().get_cds(contig_id, snp)
                    self.print_gene_stats(gene, cds, snp)
                else:
                    self.print_location(snp_analysis_obj, contig_id, snp)

            total_SNP_mutations += len(contig_evidence.SNPlist)

        self.log.info("Genome has a total of %d SNP mutations", total_SNP_mutations)

    def print_gene_stats(self, gene, cds, snp):

        if "description" in gene.qualifiers:
            description = "+".join(gene.qualifiers["description"])
            description.replace(",", " ")
        else:
            description = "no description"
        snp_type = "CDS" if cds is not None else "Intron"
        self.log.info( "%s, [%d %d %s], %s, Offset:%d, Type:%s"
                      , gene.id, gene.location.start, gene.location.end, gene.location.strand
                      , description, (snp-gene.location.start), snp_type)

    def print_location(self, snp_analysis_obj, contig_id, snp):

        upstream_gene = snp_analysis_obj.get_gene_dictionary().get_upstream_gene(contig_id, snp)
        downstream_gene = snp_analysis_obj.get_gene_dictionary().get_downstream_gene(contig_id, snp)

        if upstream_gene is not None and downstream_gene is not None:

            self.log.info( "%s, [%d %d %s], EndOffset:%d, %s, [%d %d %s], BeginOffset:%s"
                          , downstream_gene.id, downstream_gene.location.start, downstream_gene.location.end
                          , downstream_gene.location.strand, (snp-downstream_gene.location.end)
                          , upstream_gene.id, upstream_gene.location.start, upstream_gene.location.end
                          , upstream_gene.location.strand, (upstream_gene.location.start-snp))

        elif downstream_gene is not None:

            self.log.info( "%s, [%d %d %s], EndOffset:%d"
                          , downstream_gene.id, downstream_gene.location.start, downstream_gene.location.end
                          , downstream_gene.location.strand, (snp-downstream_gene.location.end))

        elif upstream_gene is not None:

            self.log.info( "%s, [%d %d %s], BeginOffset:%s"
                          , upstream_gene.id, upstream_gene.location.start, upstream_gene.location.end
                          , upstream_gene.location.strand, (upstream_gene.location.start-snp))

        else:

            self.log.warning( "No upstream or downstream gene found for contig: %s location: %d", contig_id, snp)

    def print_snp(self, contig_id, snp_index, genome_evidence):

        contig = genome_evidence[contig_id].Contigrecord
        fixed = genome_evidence[contig_id].Contigfixedarray
        insert = genome_evidence[contig_id].Contiginsertarray

        nucleotide = contig.seq[snp_index]
        nucleotide_count = fixed[snp_index][GenomeEvidence.nucleotide_offset[nucleotide]]
        a_count = fixed[snp_index][GenomeEvidence.nucleotide_offset["A"]]
        c_count = fixed[snp_index][GenomeEvidence.nucleotide_offset["C"]]
        g_count = fixed[snp_index][GenomeEvidence.nucleotide_offset["G"]]
        t_count = fixed[snp_index][GenomeEvidence.nucleotide_offset["T"]]
        delete_count = fixed[snp_index][GenomeEvidence.nucleotide_offset["-"]]
        insert_count = fixed[snp_index][GenomeEvidence.nucleotide_offset["+"]]
        count_list = [a_count, c_count, g_count, t_count, delete_count, insert_count]
        mutation_nucleotide = GenomeEvidence.nucleotide_list[count_list.index(max(count_list))]
        sum_count_list = float(sum(count_list))
        if sum_count_list > 0:
            proportion_mutant = 1.0 - (nucleotide_count / sum_count_list)
        else:
            proportion_mutant = 0.0

        self.log.info("Region:%s, Ref.Idx.Mut:%s.%d.%s, Proportion:%f, 'A':%d, 'C': %d, 'G':, %d, 'T': %d, '-':%d, '+':%d"
                     , contig.id, nucleotide, snp_index, mutation_nucleotide, proportion_mutant
                     , a_count, c_count, g_count, t_count, delete_count, insert_count)

        if max(count_list) == insert_count:
            lookback = 100
            for idx in range(snp_index-lookback, snp_index+1):
                if insert[idx] is not None:
                    for sequence, count in insert[idx].items():
                        self.log.info("+Sequence:%s, inserted:%d times, at contig location:%d", sequence, count, idx)

    def old_stuff(self, snp_index, gene, mutant_seq, contig):

        if "description" in gene.qualifiers:
            description = "+".join(gene.qualifiers["description"])
        else:
            description = "no description"

        mutant_protein = gene.extract(mutant_seq).translate()
        original_protein = gene.extract(contig.seq).translate()
        self.log.info("SNP, %d, Region, %s, Gene, %s, Function, %s", snp_index, contig.id, gene.id, description)
        if mutant_protein != original_protein:
            difference_string = ""
            for i in range(len(mutant_protein)):
                if i < len(original_protein):
                    if mutant_protein[i] == original_protein[i]:
                        difference_string += mutant_protein[i].lower()
                    else:
                        self.log.info("Mutant Residue: %s", original_protein[i]+str(i+1)+mutant_protein[i].upper())
                        difference_string += "{}=>{}".format(original_protein[i], mutant_protein[i].upper())
                else:
                    difference_string += mutant_protein[i].upper()
            self.log.info("Mutant Protein Residue Sequence  :%s", difference_string)
        else:
            self.log.info("No Protein Missense, Original Protein: %s", original_protein)

    def mutant_dna_sequence(self, gene, mutant_seq, count_list, idx):

        mutant_before = mutant_seq[:idx]
        mutant_after = mutant_seq[idx + 1:]

        if count_list[GenomeEvidence.nucleotide_offset["A"]] == max(count_list):
            mod_mutant = mutant_before + "A" + mutant_after
        elif count_list[GenomeEvidence.nucleotide_offset["C"]] == max(count_list):
            mod_mutant = mutant_before + "C" + mutant_after
        elif count_list[GenomeEvidence.nucleotide_offset["G"]] == max(count_list):
            mod_mutant = mutant_before + "G" + mutant_after
        elif count_list[GenomeEvidence.nucleotide_offset["T"]] == max(count_list):
            mod_mutant = mutant_before + "T" + mutant_after
        elif count_list[GenomeEvidence.nucleotide_offset["-"]] == max(count_list):
            mod_mutant = mutant_before + mutant_after
        else:
            self.log.warning("Unprocessed mutant, gene:%s, index %d", gene.id, idx)
            mod_mutant = mutant_seq

        return mod_mutant
