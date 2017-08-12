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
        self.gene_evidence.read_mp_sam_file(sam_file_name)
        self.count_threshold = 20

    def print_contig_stats(self):

        total_SNP_mutations = 0
        for contig_id, contig_record in self.gene_dictionary.contig_dict.items():
            contig_evidence = self.gene_evidence.get_contig_evidence(contig_id)
            contig_SNP_mutations = 0

            for gene in contig_record.Genelist:
                contig_SNP_mutations += self.print_gene_stats(gene, contig_evidence)

            self.log.info("********** Contig, %s, SNP, %d, *********", contig_id, contig_SNP_mutations)
            total_SNP_mutations += contig_SNP_mutations

        self.log.info("Genome has a total of %d SNP mutations", total_SNP_mutations)

    def print_gene_stats(self, gene, contig_evidence):

        fixed = contig_evidence.Contigfixedarray
        insert = contig_evidence.Contiginsertarray
        contig = contig_evidence.Contigrecord
        mutant_seq = contig.seq
        mutant_seq_lower = mutant_seq.lower()
        gene_SNP_mutations = 0

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
                sum_count_list = float(sum(count_list))
                if sum_count_list > 0:
                    proportion_mutant = 1.0 - (nucleotide_count / sum_count_list)
                else:
                    proportion_mutant = 0.0

                if proportion_mutant >= self.args.minMutantProportion and sum_count_list >= self.args.minMutantCount:

                    gene_SNP_mutations += 1
                    mutant_seq = self.mutant_dna_sequence(gene, mutant_seq, count_list, idx)
                    mutant_seq_lower = self.mutant_dna_sequence(gene, mutant_seq_lower, count_list, idx)
                    if False:
                        self.log.info("*** Region, %s, Gene, %s, %s, Idx, %d,"
                                      , contig.id, gene.id, description, idx)
                        self.log.info("*** Ref, %s, 'A', %d, 'C', %d, 'G', %d, 'T', %d, '-', %d, insert, %d"
                                      , nucleotide, A_count, C_count, G_count, T_count, Delete_count, Insert_count)

        if gene_SNP_mutations > 0 and False:
            mutant_protein = gene.extract(mutant_seq).translate()
            original_protein = gene.extract(contig.seq).translate()
            self.log.info("SNP, %d, Region, %s, Gene, %s, Function, %s"
                          , gene_SNP_mutations, contig.id, gene.id, description)
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

        return gene_SNP_mutations


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
