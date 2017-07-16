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
# ==============================================================================#
# MinorityReport.py
# 20150826 Jeremy Horst
# 20161220 added GTF/GFF2 option
# 20161128 added SoftClip option in CIGAR string
# 20170111 fixed CIGAR handling for multiple deletions and insertions in a read
# 20170620 modified for use by the OSM library by James McCulloch.
#
# INPUT:	FASTA of reference sequence
#			GFF gene annotation file (SAME CHROMOSOME NAMES AS FASTA!!!)
#			parent strain SAM alignment file
#			child/mutant strain SAM alignment file
#
# nonsynonymous variants
# PROCESS: *populate hash table for each position from SAM files, with evidence = counts
#           handle deletions as #-'s. 
#           handle insertions as sequence position = A[CTGCTGCTG], where A is the reference sequence position left of the insertion
#           e.g. sequence[2456821] = {'A':3,'---':426,'C':5,'T':3,'G':51, 'GTCGTACGTAGCTAGC':20}
#           check for variants over % threshold
#           check if variant is in protein coding region
#           check if variant changes the amino acid
#           compare parent to child
# OUTPUT:   report all variants over n% abundant in reads mapped to this position, unique in child w.r.t. parent
#
# copy number variants
# PROCESS:  build gene model & data structure from FASTA & GFF
#          *read in each SAM file, with the start & end of each read
#           examine tiles / sliding windows (depending on settings) across each chromosome
#          *check overlap of reads to each tile
#           take ratio of mutant to parent, in context of overall counts  
#			do statistics on distribution for read ratio for each tile
#			* = slow steps / bottlenecks
# OUTPUT:	report statistics of read ratios for each tile
# ==============================================================================#

from __future__ import print_function,  division

from BCBio.GFF.GFFParser import GFFExaminer, parse
from Bio import SeqIO, Seq
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import FeatureLocation, SeqFeature

import pprint

from OSMGFFParser import ParseGFFSeq, GffAdapter, read_gff

from math import pi, sqrt, exp, log
from collections import namedtuple
import sys
import re
import os

#############
# variables #
#############

minimum_variant_proportion = 0.9
minimum_variant_counts = 20
maximum_wildtype_proportion = 0.1
maximum_wildtype_variant_counts = 20
minimum_wildtype_total_counts = 10
desired_gene_type_in_GFF = 'cds'
GFF_description_line = 'gene'
description_key = 'description'
GTF = False
trust_nonmatching_alignment = True
position_read_report = False
verbose = False
debug = False
missense_only = True
writer = False
# cnv_p_threshold = 0.001
# cnv_ratio_threshold = 1.5
cnv = False
cnv_writer = False
window_increment = 3000
window_size = 3000
cnv_report_frequency = 3000
tile_position_mean = True
# softclips = True

gencode = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T',
           'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AAC': 'N', 'AAT': 'N', 'AAA': 'K',
           'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R', 'CTA': 'L',
           'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P',
           'CCT': 'P', 'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R',
           'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V',
           'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GAC': 'D',
           'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G',
           'GGT': 'G', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F',
           'TTT': 'F', 'TTA': 'L', 'TTG': 'L', 'TAC': 'Y', 'TAT': 'Y', 'TAA': '*',
           'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}

#####################
# general functions #
#####################
def handle_cigar(cigar, aligned_read_position):
    allowable_cigar_letters = ['MIDNSHP=X']
    # M = match or mismatch. simple; keep.
    # I = insertion. add to indexing; keep
    # D = deletion. do not count in indexing; keep.
    # N = deletion. do not count in indexing; keep.
    # S = soft clip. mismatch. Like M, but affects reported left-most alignment position when cigar stars with S; keep.
    # H = hard clip. mismatch, removes string from sequence; ignore.
    # P = padding. refers to insertion in MSA not reference; ignore.
    # = = match, same as M; keep.
    # X = mismatch, same as M; keep.
    # refseq_chromosome_ID starts AFTER softclip;
    #  find when there are no letters in front of S
    if not sum([[0, 1][h] for h in [c.isalpha() for c in cigar.split('S')[0]]]):
        aligned_read_position -= int(cigar.split('S')[0])
    insertions = []
    deletions = []
    last = 0
    read_position = 0
    for l in range(len(cigar)):
        if cigar[l] in allowable_cigar_letters:
            count = int(cigar[last:l])
            last = l
            if cigar[l] in 'MS=X':
                read_position += count
            elif cigar[l] in 'HP':
                ignore = True
            elif cigar[l] in 'I':
                insertions += [[count, read_position]]
                read_position += count
            elif cigar[l] in 'DN':
                deletions += [[count, read_position]]
    return aligned_read_position, insertions, deletions


def get_sequence_from_fasta(fasta, name):
    # return ''.join('>'.join(open(fasta).read().split('>')[1:]).split(name)[1].split('\n>')[0].split('\n')[1:])
    real_entry = '\n'  # gives blank result if real sequence not found
    for entry in open(fasta).read()[1:].split('\n>'):
        if entry.split('\n')[0].split('|')[0].strip() == name:
            real_entry = entry
            break
    return ''.join(real_entry.split('\n')[1:])




def Guassian(z):
    # approximates probability given z-score according to Gaussian
    z = float(z)
    summ = z
    value = z
    for i in range(1, 1000):
        value = value * z * z / (2 * i + 1)
        summ += value
    G = abs(0.5 + (summ / sqrt(2 * pi)) * exp(-1 * z * z / 2))
    if G > 1: G = 1 - G  # handles instances where G is barely above 1.
    return G


def rational_approximation(t):
    # Abramowitz and Stegun formula 26.2.23.
    # The absolute value of the error should be less than 4.5 e-4.
    c = [2.515517, 0.802853, 0.010328]
    d = [1.432788, 0.189269, 0.001308]
    numerator = (c[2] * t + c[1]) * t + c[0]
    denominator = ((d[2] * t + d[1]) * t + d[0]) * t + 1.0
    return t - numerator / denominator


def Guassian_inverse(p):
    assert p > 0.0 and p < 1
    if p < 0.5:
        return -1 * rational_approximation(sqrt(-2.0 * log(p)))
    else:
        return rational_approximation(sqrt(-2.0 * log(1.0 - p)))


def is_negative_strand(flag):
    # converts decimal format of flag to binary,
    # checks whether read mapped to reverse compliment.
    # i.e. determines positive vs negative strand.
    try:
        bit16 = int(bin(1)[2:][-5])
    except:
        bit16 = 0
    return bit16


def mean(numbers):
    return sum(numbers) / float(len(numbers))


def stdev(numbers, mean):
    variance = 0.0
    for i in numbers:
        variance += (mean - i) ** 2
    return sqrt(variance / len(numbers))


def z_score(x, mean, stdev):
    return (x - mean) / float(stdev)


def median(numbers):
    numbers.sort()
    return numbers[int(round(len(numbers) / 2.0))]


# ==============================================================================#


######################
# specific functions #
######################

def get_CNV_ratios(mutant_tile_reads, wt_tile_reads):
    tile_read_ratios = {}
    for chromosome in chromosomes:
        tile_read_ratios[chromosome] = {}
        len_chromo = reference_sequences[chromosome][0]
        for tile in wt_tile_reads[chromosome]:
            wt_tile_count = float(wt_tile_reads[chromosome][tile])
            mutant_tile_count = mutant_tile_reads[chromosome][tile]
            if mutant_tile_count and wt_tile_count:
                ratio = mutant_tile_count / wt_tile_count
            elif mutant_tile_count and not wt_tile_count:
                ratio = float(mutant_tile_reads[chromosome][tile]) + 1
            elif not mutant_tile_count and wt_tile_count:
                ratio = 1 / (1 + wt_tile_count)
            elif not mutant_tile_count and not wt_tile_count:
                ratio = 1 / wt_mutant_ratio
            tile_read_ratios[chromosome][tile] = log(ratio * wt_mutant_ratio, 2)
    return tile_read_ratios


def get_gene_names(tile_start, tile_end):
    # report no genes if tile occurs before first gene
    if tile_end < gff_model[chromosome][0][0]:
        gene_names = ''
    # report no genes if tile occurs after last gene
    elif tile_start >= gff_model[chromosome][-1][1]:
        gene_names = ''
    else:
        first_gene = 0
        gene_start, gene_end = gff_model[chromosome][first_gene][:2]
        # get first gene index
        while first_gene < len_gff_model_chromo - 1 and gene_end < tile_start:
            first_gene += 1
            gene_start, gene_end = gff_model[chromosome][first_gene][:2]
        start = first_gene

        # get second gene index
        last_gene = first_gene + 0
        gene_start, gene_end = gff_model[chromosome][last_gene][:2]
        while last_gene < len_gff_model_chromo - 1 and gene_end < tile_end:
            last_gene += 1
            try:
                gene_start, gene_end = gff_model[chromosome][last_gene][:2]
            except:
                if verbose: print(gene_start, gene_end, gff_model[chromosome], len(gff_model[chromosome]),
                                  file=sys.stderr)

        if gene_start <= tile_end:
            end = last_gene
        else:
            end = last_gene - 1

        gene_names = ''
        # if start > end: nada
        if start == end and start != 0:
            gene_names = gff_model[chromosome][start][9] + ' (' + gff_model[chromosome][start][3] + ')'
        elif start < end:
            for j in range(start, end + 1):
                if start != 0:
                    gene_names += gff_model[chromosome][j][9] + ' (' + gff_model[chromosome][j][3] + '); '
            gene_names = gene_names[:-2]
    return gene_names


def get_tile_counts(sequence_evidence):
    tile_reads = {}
    for chromosome in sequence_evidence:
        tile_reads[chromosome] = {}
        for tile_start in range(0, reference_sequences[chromosome][0] - window_size, window_increment):
            tile_end = tile_start + window_size
            if tile_position_mean:
                tile_reads[chromosome][tile_start] = mean(
                    [sum(sequence_evidence[chromosome][i].values()) for i in range(tile_start, tile_end)])
            else:
                tile_reads[chromosome][tile_start] = median(
                    [sum(sequence_evidence[chromosome][i].values()) for i in range(tile_start, tile_end)])
    return tile_reads


def get_counts(reads):
    # fastest of all algorithms I tested to map reads to tiles/windows of any length
    tile_reads = {}
    for chromosome in reads:
        tile_reads[chromosome] = {}
        reads_in_chromosome = len(reads[chromosome])
        if verbose:
            print('mean read length', mean([index[1] - index[0] for index in reads[chromosome]]), file=sys.stderr)
            print('chromosome', chromosome, 'reads_in_chromosome', reads_in_chromosome, file=sys.stderr)
            print('len range', len(range(0, reference_sequences[chromosome][0] - window_size, window_increment)),
                  file=sys.stderr)
        index = 0
        for tile_start in range(0, reference_sequences[chromosome][0] - window_size, window_increment):
            tile_end = tile_start + window_size

            # scoot back reads until we find a read that ends before the tile starts
            while reads[chromosome][index][1] > tile_start and index > 0:
                index -= 1

            # scoot forward reads until we find a read that ends after the tile starts
            if index < reads_in_chromosome:
                while reads[chromosome][index][1] < tile_start and index < reads_in_chromosome - 1:
                    index += 1

            # now start counting reads that overlap tile
            count = 0
            if index < reads_in_chromosome:
                while reads[chromosome][index][0] < tile_end and index < reads_in_chromosome - 1:
                    count += 1
                    index += 1

            tile_reads[chromosome][tile_start] = count
        if verbose:
            print('sum tile_reads[chromosome]', sum(tile_reads[chromosome].values()), file=sys.stderr)
    return tile_reads


def get_counts_next(reads):
    # kept for posterity. This is slower, but it is suggested to be faster throughout the interwebz.
    tile_reads = {}
    for chromosome in reads:
        tile_reads[chromosome] = {}
        reads_in_chromosome = len(reads[chromosome])
        count = 0
        for tile_start in range(0, reference_sequences[chromosome][0] - window_size, window_increment):
            tile_end = tile_start + window_size
            for read_start, read_end in reads[chromosome]:
                tile_range = set(range(tile_start, tile_end))
                read_range = set(range(read_start, read_end))
                try:
                    next(x for x in tile_range if x in read_range)
                    count += 1
                except StopIteration:
                    n = 0  # ignore failure, reward success
            tile_reads[chromosome][tile_start] = count
    return tile_reads


def does_mutation_change_amino_acid(genomic_position, variant_nucleotide, gene_model):
    gene_start, gene_end, strand, gene_id, gene_type, protein_length, AA_indices, CDS_exons, stop, description = gene_model

    # account for multi-exon genes.
    # find non-ATG starts, stop at stop codons

    # gff_model[refseq_chromosome_ID] += [[ gene_start, gene_end, strand, gene_id, gene_type, protein_length, AA_indices, CDS_exons, stop, description ]]

    # CDS_exons += [[ CDS_start, CDS_end, strand, sequence ]]

    # stop += [[ CDS_exon_i, codons.index(stop_codon) {,spillover exon length (e.g. 2 nucleotides)} ]][0]

    # AA_indices += [ [start, len(codons)], [ 3*(last_seq_offset/3.0 - last_seq_offset/3), len(codons) ], etc.]

    # negative strand: exons in reverse order, 5'-3' order stays, sequence is reverse compliment, frame is from the 3'

    # find position of variant within gene
    # find which coding region (exon) this maps to

    variant_codon = '-';
    variant_AA = '-';
    missense = '-';
    position_in_codon = 0;
    reference_codon = '-'
    this_variants_exon = False
    spillover_codon_length = False
    protein_position = 0
    gene_position = 0
    if len(stop) == 2:
        stop_exon, codon_stop = stop
    else:
        stop_exon, codon_stop, spillover_codon_length = stop
    for i in range(len(CDS_exons)):
        CDS_start, CDS_end, strand_i, sequence, frame = CDS_exons[i]
        if strand != strand_i:
            if verbose:  print('strand mismatch!', CDS_start, CDS_end, strand_i, sequence, frame, strand,
                               file=sys.stderr)
        coding_start, number_of_codons = AA_indices[i]

        # check that the last nucleotide should or should not be included (< vs <=)
        if genomic_position >= CDS_start and genomic_position < CDS_end:
            this_variants_exon = i

            # index to codon, find, translate (check if reverse compliment!!!)
            position_in_exon = genomic_position - CDS_start

            # account for negative strand
            if strand == '-':    position_in_exon = CDS_end - 1 - CDS_start - position_in_exon

            # account for actual ATG start site and inter-exon codon overlap
            position_in_coding = position_in_exon - coding_start

            # check if this position occurs before the stop codon
            if i < stop_exon or position_in_coding / 3 <= codon_stop:

                codon_start = position_in_coding / 3 * 3 + coding_start
                position_in_codon = position_in_exon - codon_start
                gene_position += position_in_coding

                # look up 3 nucleotide codon by auto-rounding with integer vs float remainder of division, and accounting for start index
                reference_codon = sequence[codon_start: (position_in_coding / 3 + 1) * 3 + coding_start]

                reference_codon_2 = reference_codon

                # handle sequence from LAST exon if in spillover
                if position_in_coding < 0:
                    reference_codon_2 = CDS_exons[i - 1][3][-3 + len(reference_codon):] + reference_codon

                # check for codon spillover between exons
                elif len(reference_codon) < 3 and len(CDS_exons) > i + 1:
                    reference_codon_2 += CDS_exons[i + 1][3][: 3 - len(reference_codon)]

                # look up AA from codon
                try:
                    reference_AA = gencode[reference_codon]
                except:
                    reference_AA = '-'
                    if verbose: print('reference codon fail:', reference_codon, file=sys.stderr)

                # handle variant positioning for codons that spillover from one exon to another

                # skip for insertions, deletions, and U's/N's
                if variant_nucleotide in ('A', 'C', 'G', 'T'):

                    try:
                        # account for negative strand by taking reverse compliments
                        # sequence already is reverse_compliment_sequence
                        if strand == '-':  variant_nucleotide = reverse_compliment(variant_nucleotide)

                        # alter sequence to contain variant
                        variant_sequence = sequence[:position_in_exon] + variant_nucleotide + sequence[
                                                                                              position_in_exon + 1:]

                        # look up 3 nucleotide codon by auto-rounding with integer division
                        variant_codon = variant_sequence[(position_in_exon - coding_start) / 3 * 3 + coding_start: (((
                                                                                                                     position_in_exon - coding_start) / 3) + 1) * 3 + coding_start]

                        variant_codon_2 = variant_codon

                        if position_in_coding < 0:
                            variant_codon_2 = CDS_exons[i - 1][3][-3 + len(variant_codon):] + variant_codon

                        # check for codon spillover between exons
                        elif len(variant_codon) < 3:
                            variant_codon_2 += CDS_exons[i + 1][3][: 3 - len(variant_codon)]

                        variant_codon = variant_codon_2

                        # look up AA from codon
                        try:
                            variant_AA = gencode[variant_codon]
                        except:
                            variant_AA = '-'

                        # check if the genetic variant results in a different amino acid === missense variation
                        if reference_AA == variant_AA:
                            missense = '-'

                        # prepare to report in standard wt_AAposition_mut format (e.g. W326R)
                        else:
                            # figure out protein position
                            protein_position += ((position_in_exon - coding_start) / 3) + 1
                            missense = reference_AA + str(protein_position) + variant_AA

                    except:
                        if verbose:
                            print(
                                'translation problem: reference_codon, genomic_position, position_in_coding, coding_start, variant_nucleotide, gene_start, gene_end, strand, len(sequence), [ i[:2] for i in CDS_exons], variant_codon',
                                file=sys.stderr)
                            print('translation problem:', reference_codon, genomic_position, position_in_coding,
                                  coding_start, variant_nucleotide, gene_start, gene_end, strand, len(sequence),
                                  [i[:2] for i in CDS_exons], variant_codon, file=sys.stderr)

        elif (genomic_position < CDS_start and strand == '+') or (genomic_position > CDS_end and strand == '-'):
            break

        else:
            protein_position += number_of_codons
            gene_position += number_of_codons * 3

    if debug: print('does_mutation_change_amino_acid:', missense, position_in_codon, reference_codon, gene_position,
                    file=sys.stderr)

    return missense, position_in_codon, reference_codon, gene_position


# @def evidence_from_sam_alignment(sequence_evidence,reference_sequences,sam_file_handle,read_mapping):
def evidence_from_sam_alignment(sequence_evidence, reference_sequences, sam_file_handle):
    print("Mapping SAM alignments from", sam_file_handle, "to reference sequence...", file=sys.stderr)

    # relate SAM alignment file to reference sequence
    # SAM format: <QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> [<TAG>:<VTYPE>:<VALUE> [...]]

    # skip header lines
    sam_file = open(sam_file_handle)
    read = sam_file.readline()
    while read[0] == '@':  read = sam_file.readline()

    # go through one read (or pair) at a time, relate to reference sequence and populate nucleotide/sequence entries
    counter = 0
    while read:

        # avoid errors from mate-pair (other paired end sequence) aligning to different chromosome
        # shouldn't trust these anyway
        refseq_chromosome_ID = read.split()[2]
        mate_refseq_chromosome_ID = read.split()[6]

        if (refseq_chromosome_ID == mate_refseq_chromosome_ID) or (mate_refseq_chromosome_ID == '=') or (
        trust_nonmatching_alignment):

            # use (I)nsertion or (D)eletion in CIGAR entry (read.split()[5])
            cigar = read.split()[5]
            if cigar != '0' and refseq_chromosome_ID != '*':

                # get alignment info
                aligned_read_position = int(read.split()[3]) - 1

                aligned_read_position, insertions, deletions = handle_cigar(cigar, aligned_read_position)

                # the MASTER script will grab nonmatching SAM lines from Paired-End reads, where the 2nd read hits the target chromosome, but not the first
                # 20170106		# PROBLEM? The mate read will be picked up on its line - no double counting!!!
                #				if trust_nonmatching_alignment:
                #					if not reference_sequences.has_key(refseq_chromosome_ID) and reference_sequences.has_key(mate_refseq_chromosome_ID):
                #						refseq_chromosome_ID = mate_refseq_chromosome_ID
                reference_sequence_length = reference_sequences[refseq_chromosome_ID][0]

                read_sequence = read.split()[9]

                read_length = len(read_sequence)

                if aligned_read_position + read_length >= reference_sequence_length:
                    aligned_read_position -= reference_sequence_length

                # iterate through query read positions
                read_index = 0
                while read_index < read_length:
                    entry = ''
                    index_adder = 0
                    variant_genomic_position = aligned_read_position + read_index

                    # insertions are added as entry into left-most nucleotide
                    for insert_length, insert_position in insertions:
                        if read_index == insert_position - 1:
                            entry = read_sequence[read_index: read_index + 1 + insert_length]
                            index_adder = insert_length

                    # deletions are added as -'s
                    for delete_length, delete_position in deletions:
                        if read_index == delete_position:
                            for deletion in range(delete_position, delete_position + delete_length):
                                entry = '-'

                                # add positions if necessary from insertions!
                                while len(sequence_evidence[refseq_chromosome_ID]) <= variant_genomic_position:
                                    sequence_evidence[refseq_chromosome_ID] += [{'A': 0, 'C': 0, 'G': 0, 'T': 0}]

                                # add gaps (-) to evidence w.r.t reference sequence / genome
                                if not sequence_evidence[refseq_chromosome_ID][variant_genomic_position].has_key(entry):
                                    sequence_evidence[refseq_chromosome_ID][variant_genomic_position][entry] = 0
                                sequence_evidence[refseq_chromosome_ID][variant_genomic_position][entry] += 1

                                # now move along w.r.t. reference sequence
                                aligned_read_position += 1
                                variant_genomic_position = aligned_read_position + read_index

                            # progress to add info for nucleotide AFTER deletion,
                            # which w.r.t the current read is the read_index
                            entry = read_sequence[read_index]

                    if entry:  # in/del's must be added to data model for this position if present
                        while len(sequence_evidence[refseq_chromosome_ID]) <= variant_genomic_position:
                            sequence_evidence[refseq_chromosome_ID] += [{'A': 0, 'C': 0, 'G': 0, 'T': 0}]
                        if not sequence_evidence[refseq_chromosome_ID][variant_genomic_position].has_key(entry):
                            sequence_evidence[refseq_chromosome_ID][variant_genomic_position][entry] = 0
                    else:
                        entry = read_sequence[read_index]

                    #					print('sequence_evidence-2',type(sequence_evidence),file=sys.stderr)
                    # store data as evidence
                    while len(sequence_evidence[refseq_chromosome_ID]) <= variant_genomic_position:
                        sequence_evidence[refseq_chromosome_ID] += [{'A': 0, 'C': 0, 'G': 0, 'T': 0}]
                    if not sequence_evidence[refseq_chromosome_ID][variant_genomic_position].has_key(entry):
                        sequence_evidence[refseq_chromosome_ID][variant_genomic_position][entry] = 0
                    sequence_evidence[refseq_chromosome_ID][variant_genomic_position][entry] += 1

                    # move on to the next
                    read_index += 1 + index_adder

                    # add read to read count
                #				insertion_length = get_insertion_length(cigar)
                #				deletion_length = get_deletion_length(cigar)

                insertion_length = sum([n[0] for n in insertions])
                deletion_length = sum([d[0] for d in deletions])

                end = aligned_read_position + read_length - insertion_length + deletion_length
                # SoftClips are not factored in to alignment (thus "soft clip")
                # @				read_mapping[refseq_chromosome_ID] += [[aligned_read_position,end]]

                counter += 1
                if not counter % 100000:
                    print('\t', ' ' * (10 - len(str(counter))) + str(counter),
                          'reads incorporated into evidence model.', file=sys.stderr)

        # load next set
        read = sam_file.readline()

    # @	# sort reads by beginning of alignment
    # @	for chromosome in read_mapping:
    # @		read_mapping[chromosome].sort()

    print('\t', counter, "SAM entries mapped", file=sys.stderr)
    # @	return sequence_evidence, read_mapping, counter
    return sequence_evidence, counter

def make_evidence_hash_table(reference_sequence_file_handle):
    print("Setting up reference evidence data structure...", file=sys.stderr)
    # set up nucleotide/sequence entries
    wt_sequence_evidence = {}
    mutant_sequence_evidence = {}
    reference_sequences = {}
    chromosomes = []

    # find all sequences
    # store names as chromosomes, get length for each
    for chromosome_entry in open(reference_sequence_file_handle).read()[1:].split('\n>'):
        chromosome_name = chromosome_entry.split()[0]
        reference_sequence = ''.join(chromosome_entry.split('\n')[1:]).upper()
        reference_sequences[chromosome_name] = [len(reference_sequence), reference_sequence]
        chromosomes += [chromosome_name]
        wt_sequence_evidence[chromosome_name] = []
        mutant_sequence_evidence[chromosome_name] = []
        print('\tReference sequence', chromosome_name, 'length:', len(reference_sequence), file=sys.stderr)
        for n in range(len(reference_sequence)):
            # set up new entry
            wt_sequence_evidence[chromosome_name] += [{'A': 0, 'C': 0, 'G': 0, 'T': 0}]
            mutant_sequence_evidence[chromosome_name] += [{'A': 0, 'C': 0, 'G': 0, 'T': 0}]

    return wt_sequence_evidence, mutant_sequence_evidence, reference_sequences, chromosomes


# ==============================================================================#

class OSMGenomeComparison(object):

    def __init__(self, args, log):

        # Shallow copies of the runtime environment.
        self.log = log
        self.args = args

    def comparison(self):

        global verbose
        global writer
        global position_read_report
        global minimum_variant_counts
        global minimum_variant_proportion
        global maximum_wildtype_proportion
        global maximum_wildtype_variant_counts
        global minimum_wildtype_total_counts
        global missense_only
        global cnv
        global cnv_writer
        global window_increment
        global window_size
        global cnv_report_frequency
        global tile_position_mean

            # load & run #
            ##############

        try:
            # if 1==1:
            # load required input arguments
            reference_sequence_file_handle = self.args.fastaFile
            gff_file_handle = self.args.gffFile
            wt_sam_file_handle = self.args.parentFile
            mutant_sam_file_handle = self.args.mutantFile

            parse_gff_seq = ParseGFFSeq(self.args, self.log)
            chr12_record = parse_gff_seq.get_id("chr12")
            ParseGFFSeq.print_gff_record("*", chr12_record)
            adapted_gff_model = GffAdapter(self.args, self.log).adapt_chromosome_seqrecord(chr12_record)

            # load optional input arguments
            if '-vp' in sys.argv:  # threshold to report positions, if the reference nucleotide is not supported by >n% of data
                minimum_variant_proportion = float(sys.argv[sys.argv.index('-vp') + 1])
                if minimum_variant_proportion > 1 or minimum_variant_proportion < 0: error
            if '-wp' in sys.argv:  # threshold to report positions, if the reference nucleotide is not supported by >n% of data
                maximum_wildtype_proportion = float(sys.argv[sys.argv.index('-wp') + 1])
                if maximum_wildtype_proportion > 1 or maximum_wildtype_proportion < 0: error
            if '-vc' in sys.argv:
                minimum_variant_counts = int(sys.argv[sys.argv.index('-vc') + 1])
                if minimum_variant_counts < 0: error
            if '-wc' in sys.argv:
                maximum_wildtype_variant_counts = int(sys.argv[sys.argv.index('-wc') + 1])
                if maximum_wildtype_variant_counts < 0: error
            if '-wtc' in sys.argv:
                minimum_wildtype_total_counts = int(sys.argv[sys.argv.index('-wtc') + 1])
                if minimum_wildtype_total_counts < 0: error
            if '-o' in sys.argv:
                writer = open(sys.argv[sys.argv.index('-o') + 1], 'w')
            if '-v' in sys.argv:
                verbose = True
                print('verbose mode', file=sys.stderr)
            if '-debug' in sys.argv:
                debug = True
                verbose = True
                print('debugging mode (extra verbose)', file=sys.stderr)
            if '-all' in sys.argv:
                missense_only = False
                print('reporting all variants whether nonsynonymous or only nucleotide changes', file=sys.stderr)
            if '-gene_type' in sys.argv:
                desired_gene_type_in_GFF = sys.argv[sys.argv.index('-gene_type') + 1]
            if ('-pe' in sys.argv):
                trust_nonmatching_alignment = False
                print('trusting only paired read alignments', file=sys.stderr)
            if '-position_read_report' in sys.argv:
                position_read_report = True
                position_read_reporter = open(sys.argv[sys.argv.index('-position_read_report') + 1], 'w')
                position_read_reporter.write('\t'.join(['chromosome', 'position', 'parent_counts', 'mutant_counts']) + '\n')
            # if '-nosoft' in sys.argv:
            #	softclips=False
            if '-description_key' in sys.argv:
                description_key = sys.argv[sys.argv.index('-description_key') + 1]
            # print('description_key changed to:',description_key, file=sys.stderr)
            if '-cnv' in sys.argv:
                cnv = True
                if '-window_size' in sys.argv:
                    window_size = int(sys.argv[sys.argv.index('-window_size') + 1])
                if '-window_increment' in sys.argv:
                    window_increment = int(sys.argv[sys.argv.index('-window_increment') + 1])
                if '-report_frequency' in sys.argv:
                    cnv_report_frequency = int(sys.argv[sys.argv.index('-report_frequency') + 1])
                if '-median' in sys.argv:
                    tile_position_mean = False
                if '-o_cnv' in sys.argv:
                    cnv_writer = open(sys.argv[sys.argv.index('-o_cnv') + 1], 'w')

                if cnv_report_frequency % window_increment:
                    print('ERROR: report_frequency must be a multiple of window_increment.', file=sys.stderr)
                    print('report_frequency:', cnv_report_frequency, '\nwindow_increment', window_increment, file=sys.stderr)
                    sys.exit()
                if window_size < window_increment:
                    print('ERROR: window_size must be equal or bigger than window_increment to achieve coverage.',
                          file=sys.stderr)
                    print('window_size:', window_size, '\nwindow_increment', window_increment, file=sys.stderr)
                    sys.exit()

            ## START MAIN ##
            print('', file=sys.stderr)
            # report variables
            if verbose:
                print('Variables are set as follows:', file=sys.stderr)
                print('\tminimum variant proportion', minimum_variant_proportion, file=sys.stderr)
                print('\tmaximum wildtype proportion', maximum_wildtype_proportion, file=sys.stderr)
                print('\tminimum variant counts', minimum_variant_counts, file=sys.stderr)
                print('\tmaximum wildtype counts', maximum_wildtype_variant_counts, file=sys.stderr)

            # load input sequence
            wt_sequence_evidence, mutant_sequence_evidence, reference_sequences, chromosomes = make_evidence_hash_table(
                reference_sequence_file_handle)

            # load gff gene model corresponding to reference sequence
            gff_model = read_gff(gff_file_handle, reference_sequences)

            adapted_list = adapted_gff_model["chr12"]
            for item in adapted_list:
                print("\nadapted:", item)
            gff_list = gff_model["chr12"]
            for item in gff_list:
                print("\ngff:", item)

            GffAdapter(self.args, self.log).compare_minority_dict(gff_model, adapted_gff_model)


            print('\t', sum([len(gff_model[chromosome]) for chromosome in chromosomes]), 'genes', file=sys.stderr)

            # load sam alignment file, align to the reference sequence
            if os.stat(mutant_sam_file_handle).st_size != 0:
                mutant_sequence_evidence, mutant_total_count = evidence_from_sam_alignment(mutant_sequence_evidence,
                                                                                           reference_sequences,
                                                                                           mutant_sam_file_handle)
                # NOTE: reference_sequence[0] corresponds to sequence_evidence[1]
                if os.stat(wt_sam_file_handle).st_size != 0:
                    wt_sequence_evidence, wt_total_count = evidence_from_sam_alignment(wt_sequence_evidence,
                                                                                       reference_sequences, wt_sam_file_handle)
                else:
                    print('ERROR: empty file:', wt_sam_file_handle, file=sys.stderr)
            else:
                print('ERROR: empty file:', mutant_sam_file_handle, file=sys.stderr)

            if verbose:
                print('SAM mapping complete.', file=sys.stderr)


        except IOError:
            sys.exit()

            # @#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
        if wt_sequence_evidence:

            # output title line
            print("Preparing output...\n", file=sys.stderr)
            printlist = ['#score', 'chromosome', 'genomic_mutation', 'gene_mutation', 'protein_mutation', 'mutation_counts',
                         'mutation_proportion', 'mutant_position_counts', '|', 'parent_mutation_counts',
                         'parent_mutation_proportion', 'parent_position_counts', '|', 'gene_id', 'strand', 'protein_length',
                         'gene_description', '||', 'A-gene_variant', 'protein_variant', 'variant_counts', 'variant_proportion',
                         '|', 'parent_variant_counts', 'parent_variant_proportion', '||', 'C-gene_variant', 'protein_variant',
                         'variant_counts', 'variant_proportion', '|', 'parent_variant_counts', 'parent_variant_proportion',
                         '||', 'G-gene_variant', 'protein_variant', 'variant_counts', 'variant_proportion', '|',
                         'parent_variant_counts', 'parent_variant_proportion', '||', 'T-gene_variant', 'protein_variant',
                         'variant_counts', 'variant_proportion', '|', 'parent_variant_counts', 'parent_variant_proportion',
                         '||', 'nonACGT-gene_variant', 'protein_variant', 'variant_counts', 'variant_proportion', '|',
                         'parent_variant_counts', 'parent_variant_proportion', '||']
            if writer:
                writer.write('\t'.join(printlist) + '\n')
            else:
                print('\t'.join(printlist), file=sys.stdout)

            # get statistics for later reporting
            mean_mutant_counts = mean(
                [sum(genomic_position.values()) for chromosome in mutant_sequence_evidence for genomic_position in
                 mutant_sequence_evidence[chromosome]])
            mean_wt_counts = mean([sum(genomic_position.values()) for chromosome in wt_sequence_evidence for genomic_position in
                                   wt_sequence_evidence[chromosome]])

            # iterate through chromosomes and positions in child (variant)
            # filter for above thresholds, compare to parent (wildtype)

            position_unique_read_counts = {}
            last_entered_position = -1
            for chromosome in chromosomes:
                gene = 0

                for genomic_position in range(reference_sequences[chromosome][0]):

                    # get reference nucleotide
                    reference_nucleotide = reference_sequences[chromosome][1][genomic_position]
                    # e.g. sequence_evidence[chromosome][2456821] = {'A':3,'---':426,'C':5,'T':3,'G':51, 'GTCGTACGTAGCTAGC':20}

                    # get sum of counts to this position
                    position_counts = float(sum(mutant_sequence_evidence[chromosome][genomic_position].values()))

                    if position_counts:

                        if position_read_report:
                            wt_position_counts = int(sum(wt_sequence_evidence[chromosome][genomic_position].values()))
                            position_read_reporter.write('\t'.join([chromosome, str(genomic_position), str(wt_position_counts),
                                                                    str(int(round(position_counts)))]) + '\n')
                        # check for variants over % threshold
                        for entry in mutant_sequence_evidence[chromosome][genomic_position]:

                            # focus on novel variants, and only report position once
                            if not entry == reference_nucleotide and not last_entered_position == genomic_position:

                                # check if there are the minimum quantity of reads
                                entry_counts = mutant_sequence_evidence[chromosome][genomic_position][entry]
                                if entry_counts >= minimum_variant_counts:

                                    # check for variants over % threshold
                                    variant_proportion = entry_counts / position_counts
                                    if variant_proportion >= minimum_variant_proportion:

                                        # wildtype data
                                        wt_position_counts = float(
                                            sum(wt_sequence_evidence[chromosome][genomic_position].values()))
                                        if wt_position_counts >= minimum_wildtype_total_counts:
                                            if wt_sequence_evidence[chromosome][genomic_position].has_key(entry):
                                                wt_entry_counts = wt_sequence_evidence[chromosome][genomic_position][entry]
                                            else:
                                                wt_entry_counts = 0
                                            wt_nt_proportion = 0
                                            if wt_position_counts:
                                                wt_nt_proportion = float(wt_entry_counts) / wt_position_counts

                                            # wt/parent read count threshold
                                            if (
                                                wt_entry_counts <= maximum_wildtype_variant_counts) and wt_nt_proportion <= maximum_wildtype_proportion:

                                                # get gene information by relating to GFF model
                                                gene_start = 0;
                                                gene_end = 0;
                                                variant_in_gene = False
                                                while gene < len(gff_model[
                                                                     chromosome]) and gene_start < genomic_position and gene_end < genomic_position:
                                                    gene_start, gene_end = gff_model[chromosome][gene][:2]
                                                    if genomic_position + 1 > gene_start and genomic_position < gene_end:
                                                        variant_in_gene = gene
                                                    gene += 1
                                                gene -= 1  # go back 1, not all, so the gene search for the next entry is quick
                                                if gene < 0: gene = 0

                                                # is the variant in the type of gene of interest (e.g. protein-encoding)?
                                                if variant_in_gene:

                                                    # get all gene details
                                                    gene_start, gene_end, strand, gene_id, gene_type, protein_length, AA_indices, CDS_exons, stop, description = \
                                                    gff_model[chromosome][variant_in_gene]

                                                    # check if genomic_position is between exons
                                                    in_exon = False
                                                    for CDS_exon in CDS_exons:
                                                        if genomic_position >= min(CDS_exon[:2]) and genomic_position < max(
                                                                CDS_exon[:2]):
                                                            in_exon = True

                                                    if in_exon:
                                                        # check if variant changes an amino acid
                                                        missense, position_in_codon, reference_codon, gene_position = does_mutation_change_amino_acid(
                                                            genomic_position, entry, gff_model[chromosome][variant_in_gene])

                                                        if (missense_only and missense != '-') or not missense_only:

                                                            ##########
                                                            # REPORT #
                                                            ##########

                                                            # chromosome
                                                            printlist = [chromosome]

                                                            ####################################
                                                            # mutant genome change
                                                            printlist += [
                                                                reference_nucleotide + str(genomic_position + 1) + entry]

                                                            # mutant gene change
                                                            # handle strand specificity to get nucleotides correct w.r.t. the gene
                                                            if strand == '-':
                                                                printlist += [reverse_compliment(reference_nucleotide) + str(
                                                                    gene_position) + reverse_compliment(entry)]
                                                            else:
                                                                printlist += [reference_nucleotide + str(gene_position) + entry]

                                                            # mutant protein change
                                                            printlist += [missense]

                                                            # entry counts, %, total
                                                            printlist += [str(entry_counts), str(round(variant_proportion, 3)),
                                                                          str(int(position_counts)), '|']

                                                            # wt/parent counts, %, total
                                                            if wt_position_counts:
                                                                printlist += [str(wt_entry_counts), str(
                                                                    round(wt_entry_counts / float(wt_position_counts), 3)),
                                                                              str(wt_position_counts), '|']
                                                            else:
                                                                printlist += [str(wt_entry_counts), 'n/a',
                                                                              str(wt_position_counts), '|']

                                                            # gene name information
                                                            printlist += [gene_id, strand, str(protein_length), description,
                                                                          '||']

                                                            ####################################
                                                            ####################################

                                                            # check non-ACGT entries first, to catch data for ACTG's
                                                            nonACGT_printlist = []
                                                            if not entry in ['A', 'C', 'G', 'T']:

                                                                # mutant data
                                                                if strand == '-':
                                                                    nonACGT_printlist += [
                                                                        reverse_compliment(reference_nucleotide) + str(
                                                                            gene_position) + reverse_compliment(entry),
                                                                        missense, str(entry_counts),
                                                                        '(' + str(round(variant_proportion, 3)) + ')', '|']
                                                                else:
                                                                    nonACGT_printlist += [
                                                                        reference_nucleotide + str(gene_position) + entry,
                                                                        missense, str(entry_counts),
                                                                        '(' + str(round(variant_proportion, 3)) + ')', '|']

                                                                # wildtype data
                                                                if wt_sequence_evidence[chromosome][genomic_position].has_key(
                                                                        entry) and wt_position_counts:
                                                                    nonACGT_printlist += [str(wt_entry_counts), '(' + str(
                                                                        round(wt_entry_counts / wt_position_counts, 3)) + ')',
                                                                                          '||']
                                                                else:
                                                                    nonACGT_printlist += ['0', '(0.000)', '||']

                                                            for nt in ['A', 'C', 'G', 'T']:

                                                                # mutant data
                                                                nt_counts = \
                                                                mutant_sequence_evidence[chromosome][genomic_position][nt]
                                                                variant_nt_proportion = nt_counts / position_counts

                                                                # wildtype data
                                                                wt_nt_counts = \
                                                                wt_sequence_evidence[chromosome][genomic_position][nt]

                                                                if strand == '+':
                                                                    variant_codon = reference_codon[
                                                                                    :position_in_codon] + nt + reference_codon[
                                                                                                               position_in_codon + 1:]
                                                                    gene_ref_nt = reference_nucleotide
                                                                    variant_nt = nt
                                                                else:
                                                                    variant_codon = reference_codon[
                                                                                    :position_in_codon] + reverse_compliment(
                                                                        nt) + reference_codon[position_in_codon + 1:]
                                                                    gene_ref_nt = reverse_compliment(reference_nucleotide)
                                                                    variant_nt = reverse_compliment(nt)

                                                                if debug:  print('is the codon working?:', nt, reference_codon,
                                                                                 gencode[reference_codon], position_in_codon,
                                                                                 variant_codon, gencode[variant_codon],
                                                                                 file=sys.stderr)

                                                                try:
                                                                    variant_AA = gencode[variant_codon]
                                                                except:
                                                                    variant_AA = '-'

                                                                # check if the genetic variant results in a different amino acid === missense variation
                                                                if missense[0] == variant_AA:
                                                                    nt_missense = '-'
                                                                else:
                                                                    for i in range(1, len(missense)):
                                                                        if missense[i].isalpha(): break
                                                                    nt_missense = missense[:i] + variant_AA

                                                                printlist += [gene_ref_nt + str(gene_position) + variant_nt,
                                                                              nt_missense, str(nt_counts),
                                                                              '(' + str(round(variant_nt_proportion, 3)) + ')',
                                                                              '|']

                                                                # wildtype data
                                                                if wt_position_counts:
                                                                    printlist += [str(wt_nt_counts), '(' + str(
                                                                        round(wt_nt_counts / float(wt_position_counts),
                                                                              3)) + ')', '||']
                                                                else:
                                                                    printlist += ['0', '(0.000)', '||']

                                                                # add importance score
                                                                # mutation_proportion * mutation_counts * mean_wt_counts / (parent_mutation_proportion * parent_mutation_counts *mean_mutant_counts)
                                                            # score = str(round( variant_proportion * (entry_counts/mean_mutant_counts)  /  ( 1 + wt_nt_proportion * (wt_entry_counts/mean_wt_counts) ),3))
                                                            score = str(round((((
                                                                                entry_counts + position_counts) / mean_mutant_counts) + (
                                                                               wt_position_counts / mean_wt_counts) - ((
                                                                                                                       position_counts - entry_counts) / position_counts)) / (
                                                                              1 + (wt_entry_counts / mean_wt_counts)), 3))
                                                            printlist = list([score]) + printlist

                                                            # report
                                                            printlist += nonACGT_printlist
                                                            if writer:
                                                                writer.write('\t'.join(printlist) + '\n')
                                                            else:
                                                                print('\t'.join(printlist), file=sys.stdout)

                                                            if debug:
                                                                print(
                                                                    '\nmissense, position in codon, reference codon, gene position, nt missense',
                                                                    missense, position_in_codon, reference_codon, gene_position,
                                                                    nt_missense, file=sys.stderr)
                                                                print(does_mutation_change_amino_acid(genomic_position, entry,
                                                                                                      gff_model[chromosome][
                                                                                                          variant_in_gene]),
                                                                      file=sys.stderr)
                                                                print(genomic_position, entry,
                                                                      gff_model[chromosome][variant_in_gene], file=sys.stderr)
                                                                print('\n', file=sys.stderr)

                                                            last_entered_position = genomic_position

            if cnv:

                # @#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
                # CNV calculations

                # calculate tile size
                #	if cnv_ratio_threshold >= 1:  TTT = Guassian_inverse(1 - 0.5*cnv_p_threshold)
                #	else:                         TTT = Guassian_inverse(0.5*cnv_p_threshold)
                #	t_threshold = abs(Guassian_inverse(cnv_p_threshold))
                #	print('\tcnv_p_threshold:',cnv_p_threshold,'t_threshold:',t_threshold,file=sys.stderr)
                #	window_size = int(round((total_mutant_reads*cnv_ratio_threshold*cnv_ratio_threshold + total_wt_reads)*reference_sequences[chromosome][0] * TTT*TTT / ((1-cnv_ratio_threshold)*(1-cnv_ratio_threshold)*total_mutant_reads*total_wt_reads)))
                #	if(window_size%2)!=0: window_size+=1  # ease window staggering
                #	print('\ttile size:',window_size,'chromosome size:',reference_sequences[chromosome][0],'ratio:',float(window_size)/reference_sequences[chromosome][0],file=sys.stderr)

                # count reads to each tile in genome for each data set
                # @	mutant_tile_reads = get_counts(mutant_reads)
                # @	wt_tile_reads = get_counts(wt_reads)
                mutant_tile_reads = get_tile_counts(mutant_sequence_evidence)
                wt_tile_reads = get_tile_counts(wt_sequence_evidence)

                # get statistics
                mean_mutant_tile_reads = mean(
                    [tile_count for chromosome in mutant_tile_reads for tile_count in mutant_tile_reads[chromosome].values()])
                mean_wt_tile_reads = mean(
                    [tile_count for chromosome in wt_tile_reads for tile_count in wt_tile_reads[chromosome].values()])

                wt_mutant_ratio = mean_wt_tile_reads / mean_mutant_tile_reads

                if verbose:
                    print('sum mutant_tile_reads', sum([x for i in mutant_tile_reads for x in mutant_tile_reads[i].values()]),
                          'sum wt_tile_reads', sum([x for i in wt_tile_reads for x in wt_tile_reads[i].values()]),
                          file=sys.stderr)
                    print('mean_mutant_tile_reads', mean_mutant_tile_reads, 'mean_wt_tile_reads', mean_wt_tile_reads,
                          file=sys.stderr)

                # get CNV ratios for all tiles
                tile_read_ratios = get_CNV_ratios(mutant_tile_reads, wt_tile_reads)

                # group together for statistical analysis
                all_tile_read_ratios = [x for i in tile_read_ratios for x in tile_read_ratios[i].values()]

                # get mean & stdev of tile count ratios
                mean_tile_read_ratios = mean(all_tile_read_ratios)
                stdev_tile_read_ratios = stdev(all_tile_read_ratios, mean_tile_read_ratios)
                if verbose:
                    print('mean, stdev of tile count ratios', mean_tile_read_ratios, stdev_tile_read_ratios, file=sys.stderr)

                # output title line
                printline = '\t'.join(
                    ['chromosome', 'start', 'end', 'wt_tile_counts', 'mutant_tile_counts', 'log2_CNV', 'probability',
                     'gene names'])
                if not cnv_writer:
                    print(printline, file=sys.stdout)
                else:
                    cnv_writer.write(printline + '\n')

                for chromosome in chromosomes:
                    len_gff_model_chromo = len(gff_model[chromosome])
                    tiles = range(0, reference_sequences[chromosome][0] - window_size, cnv_report_frequency)
                    for tile_start in tiles:

                        tile_end = tile_start + window_size

                        log2_CNV = tile_read_ratios[chromosome][tile_start]  #
                        z = z_score(log2_CNV, mean_tile_read_ratios, stdev_tile_read_ratios)

                        # log2_CNV probability calculation
                        if z >= 0:
                            p = 2 * (1 - Guassian(z))
                        else:
                            p = 2 * Guassian(z)
                        if p > 1 and verbose:  print('p>1:', p, 'log2_CNV', log2_CNV, 'z', z, 'Guassian(z)', Guassian(z),
                                                     file=sys.stderr)

                        # get gene information by relating to GFF model
                        gene_names = get_gene_names(tile_start, tile_end)

                        # look up other data for reporting
                        mutant_tile_counts = mutant_tile_reads[chromosome][tile_start]  #
                        wt_tile_counts = wt_tile_reads[chromosome][tile_start]  #

                        # report
                        printline = '\t'.join([chromosome, str(tile_start), str(tile_end), str(round(wt_tile_counts, 1)),
                                               str(round(mutant_tile_counts, 1)), str(round(log2_CNV, 3)), str(p), gene_names])
                        if not cnv_writer:
                            print(printline, file=sys.stdout)
                        else:
                            cnv_writer.write(printline + '\n')

            if writer:  writer.close()
            if position_read_report:  position_read_reporter.close()
