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
import re
import sys


# Collect the SAM file single end reads (double ended reads not supported yet.
# Nucleotide counts and deletions go into an unsigned integer numpy area the same size as the contig. region sequence.
# Insertions are entered into a numpy array of objects initially set to None. Insertions are added at the
# appropriate sequence index as dictionaries indexed by inserted sequence with count as value.


class GenomeEvidence(object):

    def __init__(self, args, log, genome_gff):
        # Shallow copies of the runtime environment.
        self.log = log
        self.args = args
        self.genome = genome_gff
        self.genome_evidence = self.genome_evidence(genome_gff)  # dict of numpy nucleotide evidence arrays
        sam_fields = "Qname Flag Rname Pos Mapquality Cigar Rnext Pnext Tlen Sequence Quality Optflags"
        self.SamRecord = namedtuple("SamRecord", sam_fields)
        self.cigar_regex = re.compile("([0-9]+[MIDNSHP=X])")
        self.CigarItem = namedtuple("CigarItem", "Code Count")
        self.nucleotide_offset = { "A" : 0, "C" : 1, "G" : 2, "T" : 3, "-" : 4}
        self.insertion = 0
        self.deletion = 0
        self.unmapped_read = 0
        self.nucleotide_mismatch = 0
        self.mismatch_threshold = 10

    def genome_evidence(self, genome_gff):

        genome_evidence = {}

        for contig_seqrecord in genome_gff:
            genome_evidence[contig_seqrecord.id] = self.make_evidence_array(contig_seqrecord)

        return genome_evidence

    def make_evidence_array(self, contig_seqrecord):

        nucleotides = 5  # for each numpy array element allocate storage for (in order} 'A', 'C', 'G', 'T' ('U'), '-'
        seq_length = len(contig_seqrecord.seq)
        numpy_shape = (seq_length, nucleotides)
        evidence_fixed_array = np.zeros(numpy_shape, np.dtype("uint32"))  # create the fixed numpy array
        numpy_shape = (seq_length,)
        evidence_insert_array = np.empty(numpy_shape, dtype=object)  # create the fixed numpy array

        self.log.info("Contig Feature id: %s, created sequence length: %d x 5 ('A','C','G','T','-')"
                      , contig_seqrecord.id, seq_length)

        return [contig_seqrecord, evidence_fixed_array, evidence_insert_array]

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

                    self.process_sam_record(sam_record)

                    if line_counter % report_increment == 0:
                        self.log.info("Processed: %d reads; Unmapped: %d, Insert: %d; Delete: %d; Mismatch: %d"
                            , line_counter, self.unmapped_read, self.insertion, self.deletion, self.nucleotide_mismatch)

        except IOError:

            self.log.error("Problem processing SAM file: %s - check the file name and directory", sam_file_name)
            sys.exit()


    def process_sam_record(self, sam_record):

        cigar_list = self.decode_cigar(sam_record.Cigar)  # returns a list of cigar tuples (Code, Count)
        current_position = int(sam_record.Pos) - 1     # Adjust for the 1 offset convention in sam files.

        if sam_record.Rname == "*":
            self.unmapped_read += 1
            return

        if sam_record.Rname not in self.genome_evidence:
            self.unmapped_read += 1
            self.log.warning("Region: %s, not found in evidence dictionary", sam_record.Rname)
            return

        evidence_insert_array = self.genome_evidence[sam_record.Rname][2]
        evidence_fixed_array = self.genome_evidence[sam_record.Rname][1]
        reference_sequence = self.genome_evidence[sam_record.Rname][0].seq
        sam_idx = 0

        for cigar in cigar_list:

            if cigar.Code not in "SH" and cigar.Count + current_position > len(reference_sequence):
                self.log.warning("Sequence Size Exceeded at Position: %d; Region: %s, len: %d, Cigar Item: (%s,%d)"
                                 , current_position, sam_record.Rname, len(reference_sequence), cigar.Code, cigar.Count)

            if cigar.Code in "MX=":

                mismatch = 0

                for idx in range(cigar.Count):

                    sam_nucleotide = sam_record.Sequence[sam_idx + idx]
                    ref_nucleotide = reference_sequence[current_position + idx]
                    if sam_nucleotide != ref_nucleotide:
                        mismatch += 1

                    offset = self.nucleotide_offset[sam_nucleotide]
                    evidence_fixed_array[current_position + idx][offset] += 1

                if mismatch >= self.mismatch_threshold:
                    self.log.info("Mismatch: %d, Region: %s, Location: %d, Ref Seq: %s, Sam Seq: %s, Cigar Item: (%s,%d), Cigar: %s, Flag: %d, OptFlag: %s"
                                  , mismatch, sam_record.Rname, current_position
                                  , reference_sequence[current_position: (current_position + cigar.Count)]
                                  , sam_record.Sequence[sam_idx: (sam_idx + cigar.Count)]
                                  , cigar.Code, cigar.Count, sam_record.Cigar, int(sam_record.Flag)
                                  , ";".join(sam_record.Optflags))

                self.nucleotide_mismatch += mismatch
                sam_idx += cigar.Count
                current_position += cigar.Count

            elif cigar.Code == "D":

                for idx in range(cigar.Count):

                    nucleotide = "-"
                    offset = self.nucleotide_offset[nucleotide]
                    evidence_fixed_array[current_position + idx][offset] += 1

                self.deletion += 1
                current_position += cigar.Count

            elif cigar.Code == "I":

                if evidence_insert_array[current_position] is None:
                    evidence_insert_array[current_position] = {}

                insert_sequence = sam_record.Sequence[sam_idx: (sam_idx+cigar.Count)]

                if insert_sequence in evidence_insert_array[current_position]:
                    evidence_insert_array[current_position][insert_sequence] += 1
#                    self.log.info("Region: %s, Location: %d, Sequence: %s previously inserted"
#                                  , sam_record.Rname, current_position, insert_sequence)
                else:
                    evidence_insert_array[current_position][insert_sequence] = 1

                sam_idx += cigar.Count
                self.insertion += 1

            elif cigar.Code == "S":
                sam_idx += cigar.Count

            elif cigar.Code == "N":
                current_position += cigar.Count

        return

    def decode_cigar(self, cigar):

        cigar_list = []
        for cigar_item in re.findall(self.cigar_regex, cigar):
            cigar_code = cigar_item[-1]
            cigar_count = int(cigar_item[:-1])
            cigar_list.append(self.CigarItem(Code=cigar_code, Count=cigar_count))

        return cigar_list
