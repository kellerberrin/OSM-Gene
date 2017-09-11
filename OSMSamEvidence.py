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

import numpy as np
from collections import namedtuple

from OSMGeneEvidence import GenomeEvidence

import libread_sam


##############################################################################################
#
#   Base class object to hold the results of reading the SAM file
#
##############################################################################################


class SamFileEvidence(object):

    EvidenceFields = namedtuple("EvidenceFields", "Contigrecord Contigfixedarray Contiginsertarray")

    def __init__(self, args, log, sam_filename):

        self.args = args
        self.log = log
        self.sam_filename = sam_filename
        self.genome_evidence = None  # Populated by the SAM file read implementation sub-class

    ##############################################################################################
    #
    #   Public class member returns the evidence object so that we can delete the sam reader.
    #
    ##############################################################################################

    def get_evidence_object(self):
        return GenomeEvidence(self.log, self.genome_evidence, self.sam_filename)


##############################################################################################
#
#   Uses a multithreaded C++ library to read and process the SAM file
#
##############################################################################################

class SamFileCPPLibrary(SamFileEvidence):

    def __init__(self, args, log, genome_gff, io_queue_size=1000000, lock_granularity=1000, processes=1, sam_filename=""):

        super(SamFileCPPLibrary, self).__init__(args, log, sam_filename)

        self.genome_evidence = self.__mp_genome_evidence(genome_gff)

        self.contig_names, self.contig_sizes = self.__mp_contig_names(genome_gff)

        self.read_sam_file_cpp()

    def __mp_contig_names(self, genome_gff):

        contig_names = []
        contig_sizes = []

        for contig_seqrecord in genome_gff:
            contig_names.append(contig_seqrecord.id)
            contig_sizes.append(len(contig_seqrecord.seq))

        return contig_names, contig_sizes

    def __mp_genome_evidence(self, genome_gff):

        genome_evidence = {}

        for contig_seqrecord in genome_gff:
            genome_evidence[contig_seqrecord.id] = self.__make_evidence_array(contig_seqrecord)

        return genome_evidence

    def __make_evidence_array(self, contig_seqrecord):

        # for each numpy array element allocate storage for (in order} 'A', 'C', 'G', 'T' ('U'), '-', '+'
        nucleotides = len(GenomeEvidence.nucleotide_list)
        seq_length = len(contig_seqrecord.seq)
        numpy_shape = (seq_length, nucleotides)
        evidence_fixed_array = np.zeros(numpy_shape, np.dtype("uint32"))  # create the fixed numpy array
        numpy_shape = (seq_length,)
        evidence_insert_array = np.empty(numpy_shape, dtype=object)  # create the insert numpy array

        self.log.info("Contig Feature id: %s, sequence length: %d", contig_seqrecord.id, seq_length)

        self.fixed = evidence_fixed_array

        return SamFileEvidence.EvidenceFields(Contigrecord=contig_seqrecord
                                             , Contigfixedarray=evidence_fixed_array
                                             , Contiginsertarray=evidence_insert_array)

    def read_sam_file_cpp(self):

        nucleotides = len(GenomeEvidence.nucleotide_list)
        numpy = libread_sam.contig_args(self.contig_names, self.contig_sizes, nucleotides)
        print(numpy[1, 1, 1])
        print(numpy[0, 0, 100])
        print(numpy[99, 999, 999])
        print(numpy[0, 0, 0])
        print("shape:", numpy.shape, " type:", numpy.dtype)

        value = 3
        self.fixed[1000000][4] = 50
        libread_sam.numpy_args(self.fixed, value, self.contig_names[0])
        print("numpy_args value:", self.fixed[3][3], " shape:", self.fixed.shape)

        for idx in range(2):
            self.leak_test(idx)

        cpp_read_lib = libread_sam.ProcessSamFile(self.args.logFilename)  # create the C++ processing class

        for contig_id, read_data in self.genome_evidence.items():
            cpp_read_lib.registerContigNumpy(contig_id, read_data.Contigfixedarray)

        cpp_read_lib.readSamFile(self.sam_filename)  # process SAM files using C++ threads

        queue_contigs = cpp_read_lib.getQueueContigs()
        queue_offsets = cpp_read_lib.getQueueOffsets()
        queue_sequences = cpp_read_lib.getQueueSequences()

    def leak_test(self, idx):

        no_leak_lib = libread_sam.NoLeak()  # test for memory leaks
        str_vec = no_leak_lib.getVec()  # return an big string
        print(" string_vec[]:", str_vec[idx], " idx:", idx)


