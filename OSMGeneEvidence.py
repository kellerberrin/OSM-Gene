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

import pysam

import pprint
import numpy as np

from math import pi, sqrt, exp, log
from collections import namedtuple
import sys
import os


class GenomeEvidence(object):

    def __init__(self, args, log, genome_gff):
        # Shallow copies of the runtime environment.
        self.log = log
        self.args = args
        self.insert_buffer_size = 1000  # overflow for  nucleotide insertions.
        self.genome_evidence = self.genome_evidence(genome_gff)  # dict of numpy nucleotide evidence arrays
        sam_fields = "Qname Flag Rname Pos Mapq Cigar Rnext Pnext Tlen Seq Qual OptF"
        self.SamRecord = namedtuple("SamRecord", sam_fields)

    def genome_evidence(self, genome_gff):

        genome_evidence = {}

        for contig_seqrecord in genome_gff:
            genome_evidence[contig_seqrecord.id] = self.make_evidence_array(contig_seqrecord)

        return genome_evidence

    def make_evidence_array(self, contig_seqrecord):

        nucleotides = 5  # for each numpy array element allocate storage for (in order} 'A', 'C', 'G', 'T' ('U'), '-'
        seq_length = len(contig_seqrecord.seq)
        numpy_shape = (seq_length + self.insert_buffer_size, nucleotides)
        evidence_array = np.zeros(numpy_shape, np.dtype("uint32"))  # create the numpy array
        self.log.info("Contig Feature id: %s, created sequence length: %d + insert buffer %d x 5 ('A','C','G','T','-')"
                      , contig_seqrecord.id, seq_length, self.insert_buffer_size)

        return [contig_seqrecord, evidence_array]

    def read_sam_file(self, sam_file_name):

        self.log.info("Processing SAM file: %s", sam_file_name)
        line_counter = 0
        report_increment = 100000

        try:

            with open(sam_file_name, 'r') as sam_handle:

                for sam_line in sam_handle:

                    line_counter += 1

                    if sam_line[0] == '@': continue  # skip header

                    sam_fields = sam_line.split("\t")

                    if len(sam_fields) < 11:
                        self.log.error("Sam file: %s, incorrect field count: %d, line num: %d, line: %s",
                                       sam_file_name, len(sam_fields), line_counter, sam_line)
                        sys.exit()

                    req_sam_fields = sam_fields[:11]
                    opt_sam_fields = [sam_fields[11:]]
                    req_sam_fields += opt_sam_fields  # should now be 12 fields, opt may be []
                    sam_record = self.SamRecord(*req_sam_fields)

                    if line_counter % report_increment == 0:
                        self.log.info("Processed: %d, reads", line_counter)

        except IOError:

            self.log.error("Problem processing SAM file: %s - check the file name and directory", sam_file_name)
            sys.exit()


