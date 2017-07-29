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
import difflib

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

    @staticmethod
    def get_CDS_list(sub_feature_list):  # recursively descend the feature tree and create a list of CDS features.

        CDS_list = []

        if sub_feature_list is not None:

            for sub_feature in sub_feature_list:

                if sub_feature.type == "CDS":

                    CDS_list.append(sub_feature)

                CDS_list = CDS_list + GeneDictionary.get_CDS_list(sub_feature.sub_features)

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
        self.count_threshold = 20

    def print_contig_stats(self):

        for contig_id, contig_record in self.gene_dictionary.contig_dict.items():
            contig_evidence = self.gene_evidence.get_contig_evidence(contig_id)
            for gene in contig_record.Genelist:
                self.print_gene_stats(gene, contig_evidence)

    def print_gene_stats(self, gene, contig_evidence):

        fixed = contig_evidence.Contigfixedarray
        insert = contig_evidence.Contifinsertarray
        contig = contig_evidence.Contigrecord
        mutant_seq = contig.seq
        mutant_seq_lower = mutant_seq.lower()
        gene_contains_mutant = False

        if "description" in gene.qualifiers:
            description = "+".join(gene.qualifiers["description"])
        else:
            description = "no description"

        cds_list = GeneDictionary.get_CDS_list(gene.sub_features)

        for cds in cds_list:

            start = cds.location.start
            end = cds.location.end

            for idx in range(start, end):

                nucleotide = contig.seq[idx]
                nucleotide_count = fixed[idx][GenomeEvidence.nucleotide_offset[nucleotide]]
                A_count = fixed[idx][GenomeEvidence.nucleotide_offset["A"]]
                C_count = fixed[idx][GenomeEvidence.nucleotide_offset["C"]]
                G_count = fixed[idx][GenomeEvidence.nucleotide_offset["G"]]
                T_count = fixed[idx][GenomeEvidence.nucleotide_offset["T"]]
                Delete_count = fixed[idx][GenomeEvidence.nucleotide_offset["-"]]
                Insert_count = 0 if insert[idx] is None else len(insert[idx])
                count_list = [A_count, C_count, G_count, T_count, Delete_count, Insert_count]

                if max(count_list) > nucleotide_count and sum(count_list) > self.count_threshold:

                    gene_contains_mutant = True
                    mutant_seq = self.mutant_dna_sequence(gene, mutant_seq, count_list, idx)
                    mutant_seq_lower = self.mutant_dna_sequence(gene, mutant_seq_lower, count_list, idx)
                    self.log.info("*** SNP Region: %s, Gene: %s CDS: %s, Idx: %d,"
                                  , contig.id, description, cds.id, idx)
                    self.log.info("*** Ref: %s, 'A': %d, 'C': %d, 'G': %d, 'T': %d, '-': %d, insert: %d"
                                  , nucleotide, A_count, C_count, G_count, T_count, Delete_count, Insert_count)

        if gene_contains_mutant:
            self.log.info("\nMutant Nucleotide Sequence: %s\n", gene.extract(mutant_seq_lower))
            mutant_protein = gene.extract(mutant_seq).translate()
            self.log.info("\nMutant Protein Amino Sequence: %s\n", mutant_protein)
            original_protein = gene.extract(contig.seq).translate()
            print('{} => {}'.format(original_protein, mutant_protein))
            for i, s in enumerate(difflib.ndiff(original_protein, mutant_protein)):
                if s[0] == " ":
                    continue
                elif s[0] == "-":
                    print(u'Delete "{}" from position {}'.format(s[-1], i))
                elif s[0] == "+":
                    print(u'Add "{}" to position {}'.format(s[-1], i))

    def mutant_dna_sequence(self, gene, mutant_seq, count_list, idx):

        mutant_before = mutant_seq[:idx]
        mutant_after = mutant_seq[idx + 1:]
        reference = mutant_seq[idx]
        print("reference:", reference)

        if count_list[0] == max(count_list):
            mod_mutant = mutant_before + "A" + mutant_after
            print("Inserted A, After:", mod_mutant[idx-10:idx+10], " Before:", mutant_seq[idx-10:idx+10])
        elif count_list[1] == max(count_list):
            mod_mutant = mutant_before + "C" + mutant_after
            print("Inserted C, After:", mod_mutant[idx-10:idx+10], " Before:", mutant_seq[idx-10:idx+10])
        elif count_list[2] == max(count_list):
            mod_mutant = mutant_before + "G" + mutant_after
            print("Inserted G, After:", mod_mutant[idx-10:idx+10], " Before:", mutant_seq[idx-10:idx+10])
        elif count_list[3] == max(count_list):
            mod_mutant = mutant_before + "T" + mutant_after
            print("Inserted T, After:", mod_mutant[idx-10:idx+10], " Before:", mutant_seq[idx-10:idx+10])
        elif count_list[4] == max(count_list):
            mod_mutant = mutant_before + mutant_after
            print("Delete Ref, After:", mod_mutant[idx-10:idx+10], " Before:", mutant_seq[idx-10:idx+10])
        else:
            mod_mutant = mutant_seq

        start = gene.location.start
        end = gene.location.end
        print("gene:", gene)
        print("\nsequence [", start, ":", end, "]:", mod_mutant[start:end], "\n")
        print("\nextract sequence:", gene.extract(mod_mutant), "\n")

        return mod_mutant
