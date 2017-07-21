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

from BCBio.GFF.GFFParser import GFFExaminer, parse
from Bio import SeqIO, Seq
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import FeatureLocation, SeqFeature

import pprint

from math import pi, sqrt, exp, log
from collections import namedtuple
import sys
import os

### Global variables used by Minority Report.

verbose = False
GFF_description_line = 'gene'
description_key = 'description'
desired_gene_type_in_GFF = 'cds'


def reverse_compliment(sequence):
    reverse_compliment_sequence = ''
    # go backwards from the end of the sequence to the beginning
    for i in range(len(sequence) - 1, -1, -1):
        if sequence[i] == 'G':
            reverse_compliment_sequence += 'C'
        elif sequence[i] == 'C':
            reverse_compliment_sequence += 'G'
        elif sequence[i] == 'A':
            reverse_compliment_sequence += 'T'
        elif sequence[i] == 'T':
            reverse_compliment_sequence += 'A'
        elif sequence[i] == '-':
            reverse_compliment_sequence += '-'
    return reverse_compliment_sequence


def prep_gene_model(CDS_exons, gene_id):
    # prepare gene for entry into gene model
    # must consider multiple exons per gene. Can check against NT & AA sequences at end of file if available.

    print("\nCDS_exons:", CDS_exons)

    # get very beginning & very end
    gene_start = min(CDS_exons[0][:2] + CDS_exons[-1][:2])
    gene_end = max(CDS_exons[0][:2] + CDS_exons[-1][:2])

    strand = CDS_exons[0][2]  # take strand of first exon
    if strand == '-':

        # swtich order of entries & sequence if sequence is on reverse strand
        CDS_exons.reverse()

        # get reverse compliment of sequence
        for i in range(len(CDS_exons)):
            CDS_exons[i][3] = reverse_compliment(CDS_exons[i][3])

    # NOTE: IF ON NEGATIVE STRAND, CDS_exons IS REVERSE COMPLIMENT
    # but the 5' and 3' arrangement for each exon remains.


    # get start/frames from GFF entry
    start = [CDS_exons[0][4], 0]
    if verbose:
        full_sequence = ''.join([i[3] for i in CDS_exons])[start[0]:]

    # get AA index for each CDS
    # AA_indices is already entered as frame, but is deciphered in the other version of this function, so we use this variable to maintain consistency outside of this function
    AA_indices = [[i[4]] for i in CDS_exons]

    # find first stop codon
    stop_codons = ('TAG', 'TAA', 'TGA')
    stop = []
    for CDS_exon_i in range(len(CDS_exons)):

        # CDS_exons data structure:
        # CDS_exons += [[ CDS_start, CDS_end, strand, sequence, frame ]]

        # check for in-translation-frame stop codons
        # divide sequence with indices...
        indexed_sequence = CDS_exons[CDS_exon_i][3][AA_indices[CDS_exon_i][0]:]

        # gather all codons in this exon as a set
        codons = [indexed_sequence[i:i + 3] for i in range(0, len(indexed_sequence), 3)]
        for stop_codon in stop_codons:
            try:
                stop += [[CDS_exon_i, codons.index(stop_codon)]]
            except ValueError:
                all_good = True

        # keep track # of codons via index
        AA_indices[CDS_exon_i] += [len(codons)]

        # check codon if one spans exons, from the end of this sequence
        spillover = int(round((len(indexed_sequence) / 3.0 - len(indexed_sequence) / 3) * 3))

        if spillover and len(CDS_exons) > CDS_exon_i + 1:
            this_exon_spillover_seq = indexed_sequence[(len(indexed_sequence) / 3) * 3:]

            next_exon_spillover_seq = CDS_exons[CDS_exon_i + 1][3][: 3 - len(this_exon_spillover_seq)]

            spillover_codon = this_exon_spillover_seq + next_exon_spillover_seq
            if not len(spillover_codon) == 3:
                if verbose:
                    print(spillover, len(CDS_exons), CDS_exon_i + 1, this_exon_spillover_seq, '.',
                          next_exon_spillover_seq, '.', file=sys.stderr)
                    print('spillover_codon is not length 3:' + spillover_codon, file=sys.stderr)

            if spillover_codon in stop_codons:
                if verbose:  print('a stop codon has been found that spans 2 different exons', file=sys.stderr)
                stop += [[CDS_exon_i, len(indexed_sequence) / 3]]

    if not stop:
        if verbose:
            print('no stop codon found for ' + gene_id, file=sys.stderr)
            print('\tCDS_exons:', CDS_exons, file=sys.stderr)
            print('\tstop_codons:', stop_codons, file=sys.stderr)
            print('\tAA_indices:', AA_indices, file=sys.stderr)
        # print('\tcodons:',[ full_sequence[i:i+3] for i in range(0,len(full_sequence),3) ],file=sys.stderr)
        #	print( ''.join(gencode[full_sequence[i:i+3]] for i in range(0,len(full_sequence),3))]
        stop = AA_indices[-1]
        protein_length = sum([i[1] for i in AA_indices])
    else:
        stop = stop[0]

        # remove exons after stop sequence
        CDS_exons = CDS_exons[:stop[0] + 1]

        # keep track of the total protein length to be able to report missense variant, e.g. M364Y where protein length = 1252 (e.g. vs 365)
        protein_length = sum([i[1] for i in AA_indices[:stop[0]]])
        protein_length += stop[1]

    return gene_start, gene_end, strand, protein_length, AA_indices, CDS_exons, stop


def read_gff(gff_file_handle, reference_sequences):

    global desired_gene_type_in_GFF

    print("Reading GFF file", file=sys.stderr)

    gene_descriptions = {}
    mRNA_parent_gene_id = {}
    line_count = 0  # +++++++

    print("gene_type:", desired_gene_type_in_GFF)   # ++++++++++

    for line in open(gff_file_handle).readlines():

        line_count += 1
#        print("++++++ sec 1, line:", line_count)

        if not line[0] == '#':
            if len(line.split()) > 4:
                entry_type = line.split()[2].lower()
                # get descriptions from gene lines
                if entry_type == GFF_description_line:
                    gene_id = line.split('ID=')[1].split(';')[0].replace('"', '').strip()
                    gene_id = gene_id.split('-')[0]
                    description = ''
                    if description_key + '=' in line:
                        description = line.split(description_key + '=')[1].split(';')[0].replace('%28', '(').replace('%29',
                                                                                                                     ')').replace(
                            '%2C', ',').replace('%2F', '/').replace('+', ' ').replace('"', '').strip()
                    elif description_key + ' ' in line:
                        description = line.split(description_key + ' ')[1].split(';')[0].replace('%28', '(').replace('%29',
                                                                                                                     ')').replace(
                            '%2C', ',').replace('%2F', '/').replace('+', ' ').replace('"', '').strip()
                    gene_descriptions[gene_id] = description
                # print(gene_id,gene_descriptions[gene_id],file=sys.stderr)

                # get mRNA name that holds splice isoform together

                if entry_type == 'mrna':  # or (not GTF and entry_type == desired_gene_type_in_GFF):
                    mRNA_id = line.split('ID=')[1].split(';')[0].replace('"', '').strip()
                    parent_gene_id = line.split('Parent=')[1].split(';')[0].replace('"', '').strip()
                    mRNA_parent_gene_id[mRNA_id] = parent_gene_id

    # load in gene starts & ends, sequence (get reverse compliment for neg strands)
    gff_model = {}
    bad_count = 0
    old_parent_mRNA = ''
    first_one = True

    line_count = 0  # ++++++
    gene_id = None
    gene_type = None
    parent_mRNA = None

    for line in open(gff_file_handle).readlines():

        line_count += 1  # +++++++

        if not line[0] == '#':

            if len(line.split()) > 4:

                gene_type = ''
                # handle multiple exons, & :. AA indexing across...
                entry_type = line.split()[2].upper()
                if 'ID=' in line:
                    gene_id = line.split('ID=')[1].split(';')[0].replace('"', '').strip()
                if 'Name=' in line:
                    gene_type = line.split('Name=')[1].split(';')[0].upper().replace('"', '').strip()
                if 'Parent=' in line:
                    parent_mRNA = line.split('Parent=')[1].split(';')[0].replace('"', '').strip()

                # non-exon containing genomes SHOULD but may not have this entry.
                print("++++++ sec 2, line:", line_count, ", gene_id:", gene_id, ", gene_type:", gene_type, ", parent:", parent_mRNA)  # +++++++

                if gene_id.lower().startswith(
                        desired_gene_type_in_GFF.lower()) \
                        or gene_type == desired_gene_type_in_GFF \
                        or entry_type == desired_gene_type_in_GFF:  # only protein coding regions will have a CDS entry
                    if parent_mRNA != old_parent_mRNA:

                        # skip for very first entry (there is no old gene), do again at end to catch last one
                        if first_one:
                            first_one = False
                        else:
                            # enter gene into gene model
                            # do AGAIN at end for very last one!
                            CDS_exons.sort()  # sometimes CDS exons are in different orders
                            gene_start, gene_end, strand, protein_length, AA_indices, CDS_exons, stop = prep_gene_model(
                                CDS_exons, old_parent_mRNA)
                            if not gff_model.has_key(refseq_chromosome_ID):  gff_model[refseq_chromosome_ID] = [
                                [0, 1, '', '', '', 0, 0, 0, 1, '']]

                            try:
                                description = gene_descriptions[mRNA_parent_gene_id[old_parent_mRNA]]
                            except:
                                description = ""

                            # remove the "rna_" at the beginning
                            if old_parent_mRNA.startswith("rna_"):
                                simplified_mRNA_name = old_parent_mRNA[4:]
                            else:
                                simplified_mRNA_name = old_parent_mRNA

                            gff_model[refseq_chromosome_ID] += [
                                [gene_start, gene_end, strand, simplified_mRNA_name, old_gene_name, protein_length,
                                 AA_indices, CDS_exons, stop, description]]

                        # reset for new gene
                        old_parent_mRNA = parent_mRNA
                        old_gene_name = gene_type

                        # reset / set even for first entry
                        CDS_exons = []

                    CDS_start = int(line.split()[3]) - 1  # to fix [1:]
                    CDS_end = int(line.split()[4])  # to fix [1:]
                    strand = line.split()[6]
                    frame = int(line.split()[7])
                    refseq_chromosome_ID = line.split()[0]
                    sequence = reference_sequences[refseq_chromosome_ID][1][CDS_start: CDS_end]
                    CDS_exons += [[CDS_start, CDS_end, strand, sequence, frame]]

    print("++++++ For Loops Satisfied ++++++")

# enter CDS entry line into gene model AGAIN at end for very last entry!

    print("last loop", line_count, ", gene_id:", gene_id, ", gene_type:", gene_type, ", parent:", parent_mRNA)  # +++++++

    if gene_id.lower().startswith(
            desired_gene_type_in_GFF.lower()) \
            or gene_type == desired_gene_type_in_GFF \
            or entry_type == desired_gene_type_in_GFF:

        print("++++++ If Statement True ++++++")

        gene_start, gene_end, strand, protein_length, AA_indices, CDS_exons, stop = prep_gene_model(CDS_exons, parent_mRNA)

        print("++++++ If Statement Satisfied ++++++")

        if not gff_model.has_key(refseq_chromosome_ID):  gff_model[refseq_chromosome_ID] = []

        try:
            description = gene_descriptions[mRNA_parent_gene_id[parent_mRNA]]
        except:
            description = ''

        simplified_mRNA_name = parent_mRNA[4:]  # remove the "rna_" at the beginning
        gff_model[refseq_chromosome_ID] += [
            [gene_start, gene_end, strand, simplified_mRNA_name, gene_type, protein_length, AA_indices, CDS_exons, stop,
             description]]

    print("about to type gf_model")
    print("++++++ Return Called ++++++", type(gff_model), len(gff_model))

    for key, value in gff_model.iteritems(): # +++++++++
        print ("gff_model_key:", key, "value type length", type(value), len(value))

    return gff_model

######################################################################################################################
##
## Biopythin based GFF and fasta parser.
##
##
######################################################################################################################

class ParseGFFSeq(object):  # parses the input gff(3) file and annotates it with a fasta sequence

    def __init__(self, args, log):
        # Shallow copies of the runtime environment.
        self.log = log
        self.args = args
        # Parse the fasta and gff files.
        self.fasta_sequence = self.__read_fasta(self.args.fastaFile)
        self.parsed_structure = self.__gff_parser(self.args.gffFile, base_dict=self.fasta_sequence)

    def __read_fasta(self, fasta_filename):

        with open(fasta_filename, 'r') as fasta_handle:
            seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))
        for key, seq in seq_dict.items():
            seq.seq.alphabet = DNAAlphabet()  # The fasta file should contain DNA code.
            self.log.info("Read DNA Fasta ID: %s, length: %d", seq.id, len(seq.seq))
        return seq_dict


    def __gff_parser(self, gff_file_name, base_dict=None):

        with open(gff_file_name) as gff_handle:
            parsed_gff = list(parse(gff_handle, base_dict=base_dict, limit_info=None))
            for rec in parsed_gff:
                self.log.info("Read Gff ID: %s, features: %d", rec.id, len(rec.features))

        return parsed_gff

    def __gff_struct_print(self, gff_file_name):

        with open(gff_file_name) as gff_handle:
            examiner = GFFExaminer()
            pprint.pprint(examiner.available_limits(gff_handle))

    def get_id(self, ident):

        for record in self.parsed_structure:

            if ident == record.id:

                return record

        return None

    def print_gff_seq(self, ident="*"):

        print("\n****************** Fasta Structure **************************************")
        print("Fasta type:", type(self.fasta_sequence), "\nFasta record:", self.fasta_sequence)

        print("\n****************** Parsed Structure Dictionary **************************************")
        print("Parsed type:", type(self.parsed_structure), "\nParsed record:", self.parsed_structure)

        for record in self.parsed_structure:
            print("\nParsed Record type:", type(record))
            ParseGFFSeq.print_gff_record(ident, record)

    @staticmethod
    def print_gff_record(ident, record):

        print("\n****************** Printing Parsed GFF Structure **************************************")
        print("Type:", type(record), "\nRecord:", record)

        for feature in record.features:
            ParseGFFSeq.__print_feature(ident, record.seq, feature, 1)

    @staticmethod
    def __print_feature(ident, seq, feature, level):

        if ident == "*" or ident in feature.id:

            print("\n******* feature level:", level, "\n", feature)

            if feature.type == "CDS":
                print("\nNucleotide Sequence:", feature.extract(seq))
                print("\nProtein Amino Sequence:", feature.extract(seq).translate())

            for sub_feature in feature.sub_features:
                ParseGFFSeq.__print_feature(ident, seq, sub_feature, level + 1)  # recursively print sub_features.


class GffAdapter(object):

    def __init__(self, args, log):

        # Shallow copies of the runtime environment.
        self.log = log
        self.args = args

        # define some named tuples as records used to pass data around.

        self.Extended_AA_Indices = namedtuple("Extended_AA_Indices", "Phase ProteinLen NucleotideLen")
        self.CDS_exons = namedtuple("CDS_exons", "Start End Strand Sequence Phase")
        self.Feature = namedtuple("Feature", "GeneStart GeneEnd Strand Id Type ProteinLen AA_indices CDS_exons Stop Description")

    def adapt_chromosome_seqrecord(self, record):  # Adapt a SeqRecord object to a list of minority style records
        feature_record_list = []

        for feature in record.features:  # process all the features in a biopython seqrecord object.
            feature_record = self.process_feature(feature, record.seq)  # returned as a namedtuple: Feature
            feature_record_list.append(list(feature_record))  # convert to a list of lists

        return {record.id: feature_record_list}  # convert to a dictionary indexed by the chromosome name.

    def process_feature(self, feature, seq):

        if feature.location.strand is None:
            self.log.error("Undefined strand for feature id:%s", feature.id)
            sys.exit()

        strand = "+" if feature.location.strand > 0 else "-"
        id = feature.id
        type = "CDS"
        CDS_exons, AA_indices, nucleotide_length, CDS_start, CDS_end = self.process_CDS_exons(feature.sub_features, seq)

        AA_indices = self.convert_AA_indices(AA_indices)

        if nucleotide_length % 3 != 0:
            self.log.warning("Feature id: %s, nucleotide length: %d is not modulo 3", feature.id, nucleotide_length)

        stop = self.process_stop(id, nucleotide_length, strand, CDS_exons)

        protein_length = int(nucleotide_length // 3) -1

        if "description" in feature.qualifiers:
            description = ";".join(feature.qualifiers["description"])
        else:
            description = "no feature description"

        CDS_exon_list = [list(CDS) for CDS in CDS_exons]  # Convert a list of namedtuples to a list of lists

        return self.Feature(GeneStart=CDS_start, GeneEnd=CDS_end, Strand=strand, Id=id, Type=type, ProteinLen=protein_length
                            , AA_indices=AA_indices, CDS_exons=CDS_exon_list, Stop=stop, Description=description)

    def convert_AA_indices(self, AA_indices):  # List of named tuple Extended_AA_Indices

        converted_AA_indices = []
        length_AA_indices = len(AA_indices)

        for idx, CDS_index in enumerate(AA_indices):
            adjusted_nucleotide_length = CDS_index.NucleotideLen - CDS_index.Phase

            if adjusted_nucleotide_length % 3 != 0:

                if idx + 1 < length_AA_indices:

                    adjusted_nucleotide_length += AA_indices[idx+1].Phase

            protein_length = int(adjusted_nucleotide_length // 3)
            converted_AA_indices.append([CDS_index.Phase, protein_length])

        return converted_AA_indices

    def get_CDS_list(self, sub_feature_list):  # recursively descend the feature tree and create a list of CDS features.

        CDS_list = []

        if sub_feature_list is not None:

            for sub_feature in sub_feature_list:

                if sub_feature.type == "CDS":

                    CDS_list.append(sub_feature)

                CDS_list = CDS_list + self.get_CDS_list(sub_feature.sub_features)

        return CDS_list

    def process_CDS_exons(self, sub_feature_list, seq):   # recursively process all gene subfeatures looking for 'CDS'

        CDS_exons_list = []
        AA_indicies_list = []
        sum_nucleotide_length = 0
        CDS_start = 0
        CDS_end = 0

        for CDS_feature in self.get_CDS_list(sub_feature_list):

            start = int(CDS_feature.location.start)
            end = int(CDS_feature.location.end)

            if CDS_start == 0 or start < CDS_start:
                CDS_start = start

            if CDS_end == 0 or end > CDS_end:
                CDS_end = end

            if CDS_feature.location.strand is None:
                self.log.error("Undefined strand for feature id: %s", CDS_feature.id)
                sys.exit()

            strand = "+" if CDS_feature.location.strand > 0 else "-"
            sequence = str(CDS_feature.extract(seq))

            if "phase" in CDS_feature.qualifiers:

                if len(CDS_feature.qualifiers["phase"]) == 1:
                    phase = int(CDS_feature.qualifiers["phase"][0])
                else:
                    self.log.error("Feature id: %s; Unexpected size for qualifier 'phase' list: %d, expected 1"
                                   , CDS_feature.id, len(CDS_feature.qualifiers["phase"]))
                    sys.exit()

            else:

                self.log.error("Feature id: %s; Expected qualifier 'phase'", CDS_feature.id)
                sys.exit()

            if "size" in CDS_feature.qualifiers:

                if len(CDS_feature.qualifiers["size"]) == 1:

                    nucleotide_length = int(CDS_feature.qualifiers["size"][0])

                else:
                    self.log.error("Feature id: %s; Unexpected size for qualifier 'size' list: %d, expected 1"
                                   , CDS_feature.id, len(CDS_feature.qualifiers["size"]))
                    sys.exit()

            else:

                self.log.error("Feature id: %s; Expected qualifier 'size' in 'CDS' record", CDS_feature.id)
                sys.exit()

            sum_nucleotide_length += nucleotide_length
            protein_length = int((nucleotide_length-phase) // 3)
            AA_indicies_list.append(self.Extended_AA_Indices(Phase=phase, ProteinLen=protein_length
                                                             ,  NucleotideLen=nucleotide_length))
            CDS_exons_list.append(self.CDS_exons(Start=start, End=end, Strand=strand
                                                 , Sequence=sequence, Phase=phase))

            print("\n++++++++++++sequence type:", type(sequence), "sequence:", sequence, "\n\n")

        return CDS_exons_list, AA_indicies_list, sum_nucleotide_length, CDS_start, CDS_end

    def process_stop(self, id, nucleotide_length, strand, CDS_exons):   # Examine all sequences looking for a stop codon and ensure order.

        stop_codons = ("TAG", "TAA", "TGA")
        start_codon = "ATG"
        stop_idx = []
        stop_idx_1 = -1
        stop_idx_2 = -1
        stop_idx_3 = -1

        if not CDS_exons:  # feature contains no sequences (chromosome feature).
            return stop_idx

        for idx, exon in enumerate(CDS_exons):

            phase = exon.Phase
            search_seq = exon.Sequence[phase:]

            if idx >= 1:  # compare the ordering of the sequences.

                previous_end = CDS_exons[idx-1].End  # The end index of the previous sequence.
                current_start = CDS_exons[idx].Start  # The start index of the current sequence

                if strand == "+" and previous_end >= current_start:
                    self.log.error("Feature id: %s; '+' strand, previous sequence end: %d, current start: %d"
                                   , id, previous_end, current_start)
                    sys.exit()

                elif strand == "-" and previous_end <= current_start:
                    self.log.error("Feature id: %s; '-' strand, previous sequence end: %d, current start: %d"
                                   , id, previous_end, current_start)
                    sys.exit()

            else:  # is the first sequence so try and find the start codon

                start_codon_offset = search_seq.find(start_codon)

                if start_codon_offset > 0:
                    self.log.warning("Feature id: %s, phase: %d adjusted sequence: %d length: %d, start codon at Nucleotide offset %d,"
                                     , id, phase, idx, len(search_seq), start_codon_offset)

                elif start_codon_offset == -1:
                    self.log.warning("Feature id: %s, phase: %d adjusted sequence: %d length: %d, start codon Not Found"
                                     , id, phase, idx, len(search_seq))

            if len(search_seq) % 3 != 0 and nucleotide_length % 3 != 0:
                self.log.warning("Feature id: %s, phase: %d adjusted sequence: %d length: %d is not modulo 3"
                                 , id, phase, idx, len(search_seq))

            search_codons = int(len(search_seq) // 3)

            for codon_idx in range(0, search_codons):

                codon_str_idx = codon_idx * 3
                codon = search_seq[codon_str_idx: codon_str_idx + 3]

                if codon == stop_codons[0]:
                   stop_idx_1 = codon_idx

                if codon == stop_codons[1]:
                   stop_idx_2 = codon_idx

                if codon == stop_codons[2]:
                   stop_idx_3 = codon_idx

            if (stop_idx_1 >= 0 and stop_idx_2 >= 0) \
                or (stop_idx_1 >= 0 and stop_idx_3 >= 0) \
                or (stop_idx_2 >= 0 and stop_idx_3 >= 0):

                self.log.warning("Feature id: %s; multiple stop codons seq phase: %d, sequence: %d, codon locations: ['TAG' %d, 'TAA' %d, 'TGA' %d]"
                                 , id, phase, idx, stop_idx_1, stop_idx_2, stop_idx_3)

            stop_idx_0 = max(stop_idx_1, stop_idx_2, stop_idx_3)

            if stop_idx_0 >= 0 and stop_idx:
                self.log.warning("Feature id: %s; stop codon found in sequence: %d, previously found in sequence: %d"
                                 , id, idx, stop_idx[0])

            if stop_idx_0 >= 0:
                stop_idx = [idx, stop_idx_0]

        if not stop_idx:
            self.log.warning("Feature id: %s; no stop codon found, last codon on last sequence assumed", id)
            stop_idx = [len(CDS_exons)-1, int(len(CDS_exons[-1].Sequence[CDS_exons[-1].Phase:]) // 3)-1]

        if stop_idx[0] != (len(CDS_exons)-1):
            self.log.warning("Feature id: %s; stop codon not found in last sequence, found in seq index %d, codon %d"
                             , id, stop_idx[0], stop_idx[1])

        return stop_idx

    def compare_minority_records(self, feature_list_a, feature_list_b):

        compare_result = True
        field_list = ["gene_start", "gene_end", "strand", "ident_name", "gene_type", "protein_length"
                      , "AA_indices", "CDS_exons", "stop", "description"]

        if len(feature_list_a) != 10 or len(feature_list_b) != 10:
            self.log.error("Compare records, Minority records are an incorrect length")
            compare_result = False

        feature_a = self.Feature(*feature_list_a)  # convert list based feature record to namedtuple
        feature_b = self.Feature(*feature_list_b)

        if feature_a.Id not in feature_b.Id and feature_b.Id not in feature_a.Id:
            self.log.error("Compare minority ident not equal record a: %s, record b: %s", feature_a.Id, feature_b.Id)
            compare_result = False

        for idx in range(len(feature_a)):

            if idx != 3 and idx != 9:

                if feature_a[idx] != feature_b[idx]:
                    self.log.error("Record a: %s, record b: %s Inequality in field %s"
                                   , feature_a.Id, feature_b.Id, field_list[idx])
                    compare_result = False

        return compare_result

    def compare_minority_dict(self, dict_a, dict_b):  # Assumes a chromosome id dictionary with a list of list features.

        error_count = 0
        max_error_count = 200
        compare_result = True

        for key, value in dict_a.iteritems():
            if key not in dict_b:
                self.log.error("Chromosome key: %s in dict_a not in dict_b", key)
                sys.exit()

        for key, value in dict_b.iteritems():
            if key not in dict_a:
                self.log.error("Chromosome key: %s in dict_b is not in dict_a", key)
                sys.exit()

        for key, value in dict_a.iteritems():

            list_a = value  # extract a list of features from the dictionary.
            list_b = dict_b[key]
            self.log.info("Comparing chromosome: %s", key)

            if len(list_a) != len(list_b):
                self.log.error("Chromosome: %s feature lists are different size, len(list_a): %d, len(list_b): %d"
                                 , key, len(list_a), len(list_b))

                for feature_a_list in list_a:
                    feature_a = self.Feature(*feature_a_list)  # convert from a list to a namedtuple
                    not_found = True
                    for feature_b_list in list_b:
                        feature_b = self.Feature(*feature_b_list)  # convert from a list to a namedtuple
                        if len(feature_a.Id) > 0 and len(feature_b.Id) > 0\
                                and (feature_a.Id in feature_b.Id or feature_b.Id in feature_a.Id):
                            not_found = False
                            if not self.compare_minority_records(feature_a, feature_b):
                                error_count += 1
                                compare_result = False

                            if error_count >= max_error_count:
                                self.log.error("Maximum: %d comparison errors reached", max_error_count)
                                sys.exit()
                            break

                    if not_found:
                        self.log.error("Feature id: %s found in dict_a, Not Found in dict_b", feature_a[3])

                for feature_b_list in list_b:
                    feature_b = self.Feature(*feature_b_list)  # convert from a list to a namedtuple
                    not_found = True
                    for feature_a_list in list_a:
                        feature_a = self.Feature(*feature_a_list)  # convert from a list to a namedtuple
                        if len(feature_a.Id) > 0 and len(feature_b.Id) > 0\
                                and (feature_a.Id in feature_b.Id or feature_b.Id in feature_a.Id):
                            not_found = False
                            if not self.compare_minority_records(feature_a, feature_b):
                                error_count += 1
                                compare_result = False

                            if error_count >= max_error_count:
                                self.log.error("Maximum: %d comparison errors reached", max_error_count)
                                sys.exit()
                            break

                    if not_found:
                        self.log.error("Feature id: %s found in dict_b, Not Found in dict_a", feature_b[3])

        return compare_result


