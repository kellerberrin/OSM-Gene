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

from BCBio.GFF.GFFParser import parse
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet

######################################################################################################################
##
## Biopython based GFF and fasta parser.
##
######################################################################################################################

class ParseGFFSeq(object):  # parses the input gff(3) file and annotates it with a fasta sequence

    def __init__(self, log, fasta_filename, gff_file_name):
        # Shallow copies of the runtime environment.
        self.log = log
        # Parse the fasta and gff files.
        self.fasta_sequence = self.__read_fasta(fasta_filename)
        self.parsed_structure = self.__gff_parser(gff_file_name, base_dict=self.fasta_sequence)

    def get_parsed_structure(self):

        return self.parsed_structure

    def get_fasta_sequence(self):

        return self.fasta_sequence

    def get_id(self, contig_id):

        for contig in self.parsed_structure:
            if contig_id == contig.id:
                return contig
        return None

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

