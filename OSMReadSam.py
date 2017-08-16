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
import multiprocessing
from multiprocessing.sharedctypes import RawArray
from ctypes import c_uint32, c_bool
from collections import namedtuple
import re
import sys

from OSMGeneEvidence import GenomeEvidence


##############################################################################################
#
#   A private object to trade-off between memory usage and lock granularity (speed)
#
##############################################################################################


class GranularityLock(object):

    def __init__(self, array_size, granularity):

        self.array_size = array_size
        self.granularity = granularity
        self.lock_count = int(array_size / granularity) + 1
        numpy_shape = (self.lock_count,)
        self.lock_array = np.empty(numpy_shape, dtype=object)  # create the lock numpy array
        for idx in range(self.lock_count):
            self.lock_array[idx] = multiprocessing.Lock()

    def acquire(self, array_index):
        idx = int(array_index/self.granularity)
        self.lock_array[idx].acquire()

    def release(self, array_index):
        idx = int(array_index/self.granularity)
        self.lock_array[idx].release()


##############################################################################################
#
#   Processes SAM files using multiple processes for maximum speed
#
##############################################################################################


class ReadSamFile(object):

    SamRecord = namedtuple("SamRecord", "Qname Flag Rname Pos Mapquality Cigar Rnext Pnext Tlen Sequence Quality Optflags")
    EvidenceFields = namedtuple("EvidenceFields", "Contigrecord Contigfixedarray Contiginsertarray Contiginsertlist Contiglocks")
    InsertTuple = namedtuple("InsertTuple", "Contigid Sequenceidx Insertseq")
    CigarItem = namedtuple("CigarItem", "Code Count")

    def __init__(self, log, genome_gff, io_queue_size=1000000, lock_granularity=1000, processes=1, sam_filename=""):
        # Shallow copies of the runtime environment.
        self.log = log
        self.genome = genome_gff
        self.lock_granularity = lock_granularity
        self.io_queue_size = io_queue_size
        # Constants and named tuple definitions
        self.cigar_regex = re.compile("([0-9]+[MIDNSHP=X])")
        self.eof_flag = -1
        self.report_increment = 100000
        # The variables below are shared between processes
        self.genome_evidence = self.__genome_evidence(genome_gff)  # dict of numpy nucleotide evidence arrays
        self.line_counter = multiprocessing.Value('i', 0)
        self.insertion = multiprocessing.Value('i', 0)
        self.deletion = multiprocessing.Value('i', 0)
        self.unmapped_read = multiprocessing.Value('i', 0)
        self.nucleotide_mismatch = multiprocessing.Value('i', 0)
        self.eof = multiprocessing.Value(c_bool, False)
        # Read the SAM file
        self.__evidence_from_sam_file(sam_filename, processes)

    ##############################################################################################
    #
    #   Public class members
    #
    ##############################################################################################

    def get_evidence_object(self):
        return GenomeEvidence(self.log, self.genome_evidence)

    ##############################################################################################
    #
    #   Private class members
    #
    ##############################################################################################


    def __evidence_from_sam_file(self, sam_filename, processes):

        self.__read_mp_sam_file(sam_filename, processes)

        self.log.info("****Completed Processing: %d reads; Unmapped: %d, Insert: %d; Delete: %d; Mismatch: %d"
                      , self.line_counter.value, self.unmapped_read.value, self.insertion.value
                      , self.deletion.value, self.nucleotide_mismatch.value)

    def __read_mp_sam_file(self, sam_file_name, child_processes):

        manager = multiprocessing.Manager()
        sam_record_queue = multiprocessing.Queue(self.io_queue_size)
        insert_record_queue = manager.Queue()
        counter_lock = multiprocessing.Lock()
        # start the SAM analysis processes.
        consumer_processes = []
        for _ in range(child_processes): # start the SAM record consumer processes
            process = multiprocessing.Process(target=self.__consume_sam_records
                                              , args=(sam_record_queue, insert_record_queue, counter_lock))
            process.daemon = True
            process.start()  # process SAM records.
            consumer_processes.append(process)

        self.log.info("Started: %d SAM record analysis processes", child_processes)

        # produce records for the consumer processes
        self.__read_sam_file(sam_file_name, sam_record_queue, child_processes)

        for idx in range(len(consumer_processes)):
            consumer_processes[idx].join()  # wait for the consumer processes to complete

        # Unpack the insert queue
        self.__process_insert_queue(insert_record_queue)

    def __genome_evidence(self, genome_gff):

        genome_evidence = {}

        for contig_seqrecord in genome_gff:
            genome_evidence[contig_seqrecord.id] = self.__make_evidence_array(contig_seqrecord)

        return genome_evidence

    def __make_evidence_array(self, contig_seqrecord):

        nucleotides = 5  # for each numpy array element allocate storage for (in order} 'A', 'C', 'G', 'T' ('U'), '-'
        seq_length = len(contig_seqrecord.seq)
        numpy_shape = (seq_length, nucleotides)
        evidence_fixed_array = self.__get_shared_numpy(numpy_shape)  # create the fixed numpy array
        numpy_shape = (seq_length,)
        evidence_insert_array = np.empty(numpy_shape, dtype=object)  # create the insert numpy array

        self.log.info("Contig Feature id: %s, sequence length: %d", contig_seqrecord.id, seq_length)

        return ReadSamFile.EvidenceFields(Contigrecord=contig_seqrecord
                                         , Contigfixedarray=evidence_fixed_array
                                         , Contiginsertarray=evidence_insert_array
                                         , Contiginsertlist=[]
                                         , Contiglocks=GranularityLock(seq_length, self.lock_granularity))

    def __get_shared_numpy(self, numpy_shape):  # The fixed evidence array is shared between processes

        shared_array = RawArray(c_uint32, (numpy_shape[0] * numpy_shape[1]))
        dt = np.dtype("uint32")
        flat_np = np.frombuffer(shared_array, dt)
        shared_np = np.reshape(flat_np, numpy_shape)

        return shared_np

    def __process_insert_queue(self, insert_record_queue):

#        insert_record_queue.put(self.eof_flag)  # insert an eof marker at the back of the queue
        insert_record_counter = 0

        while not insert_record_queue.empty():  # unpack the queue until the eof_flag is reached.
            insert_list = insert_record_queue.get()
            insert_tuple = ReadSamFile.InsertTuple( Contigid=insert_list[0]
                                                  , Sequenceidx=insert_list[1]
                                                  , Insertseq=insert_list[2])
            contig_evidence = self.genome_evidence[insert_tuple.Contigid]
            evidence_insert_array = contig_evidence.Contiginsertarray

            if evidence_insert_array[insert_tuple.Sequenceidx] is None:
                insert_dict = {}
                evidence_insert_array[insert_tuple.Sequenceidx] = insert_dict  # an indexed array of insertion dicts
            if insert_tuple.Insertseq in evidence_insert_array[insert_tuple.Sequenceidx]:
                evidence_insert_array[insert_tuple.Sequenceidx][insert_tuple.Insertseq] += 1
            else:
                evidence_insert_array[insert_tuple.Sequenceidx][insert_tuple.Insertseq] = 1

            insert_record_counter += 1

        self.log.info("Processed: %d inserted nucloetide sequences", insert_record_counter)

    def __read_sam_file(self, sam_file_name, sam_record_queue, process_count):

        self.log.info("Processing SAM file: %s", sam_file_name)

        try:

            with open(sam_file_name, 'r') as sam_handle:

                for sam_line in sam_handle:

                    sam_record_queue.put(sam_line)

            for _ in range(process_count):
                sam_record_queue.put(self.eof_flag)

        except IOError:

            for _ in range(process_count):
                sam_record_queue.put(self.eof_flag)
            self.log.error("Problem processing SAM file: %s - check the file name and directory", sam_file_name)
            sys.exit()

    def __consume_sam_records(self, sam_record_queue, insert_record_queue, counter_lock):

        while not self.eof.value:

            sam_line = sam_record_queue.get()
            if sam_line == self.eof_flag:
                self.eof.value = True
                break

            if sam_line[0] == '@': continue  # skip header

            sam_fields = sam_line.split("\t")
            if len(sam_fields) < 11:
                self.log.error("Incorrect field count: %d, line num: %d, line: %s"
                               , len(sam_fields), self.line_counter, sam_line)
                sys.exit()

            req_sam_fields = sam_fields[:11]
            opt_sam_fields = [sam_fields[11:]]
            req_sam_fields += opt_sam_fields  # should now be 12 fields, opt may be []
            sam_record = ReadSamFile.SamRecord(*req_sam_fields)

            self.__process_sam_record(sam_record)

            with counter_lock:
                self.line_counter.value += 1
                if self.line_counter.value % self.report_increment == 0:
                    self.log.info("Processed: %d reads; Unmapped: %d, Insert: %d; Delete: %d; Mismatch: %d"
                                  , self.line_counter.value, self.unmapped_read.value, self.insertion.value
                                  , self.deletion.value, self.nucleotide_mismatch.value)

        self.__queue_insert_records(insert_record_queue)   # inserted nucleotides are processed in the mainline process.

    def __process_sam_record(self, sam_record):

        cigar_list = self.__decode_cigar(sam_record.Cigar)  # returns a list of cigar tuples (Code, Count)
        current_position = int(sam_record.Pos) - 1     # Adjust for the 1 offset convention in sam files.

        if sam_record.Rname == "*":
            self.unmapped_read.value += 1
            return

        if sam_record.Rname not in self.genome_evidence:
            self.unmapped_read.value += 1
            self.log.error("Contiguous Region: %s not found in gff, probable SAM/gff mismatch", sam_record.Rname)
            sys.exit()

        contig_evidence = self.genome_evidence[sam_record.Rname]
        evidence_fixed_array = contig_evidence.Contigfixedarray
        evidence_insert_list = contig_evidence.Contiginsertlist
        reference_sequence = contig_evidence.Contigrecord.seq
        granularity_locks = contig_evidence.Contiglocks

        sam_idx = 0

        for cigar in cigar_list:

            if cigar.Code not in "SH" and cigar.Count + current_position > len(reference_sequence):
                self.log.warning("Sequence Size Exceeded at Position: %d; Region: %s, len: %d, Cigar Item: (%s,%d)"
                                 , current_position, sam_record.Rname, len(reference_sequence), cigar.Code, cigar.Count)

            if cigar.Code in "MX=":

                mismatch = 0

                for idx in range(cigar.Count):

                    sam_nucleotide = sam_record.Sequence[sam_idx + idx]
                    seq_index = current_position + idx
                    ref_nucleotide = reference_sequence[seq_index]

                    if sam_nucleotide in GenomeEvidence.nucleotide_offset:
                        offset = GenomeEvidence.nucleotide_offset[sam_nucleotide]
                        # Shared memory is protected with a granularity lock
                        granularity_locks.acquire(seq_index)
                        evidence_fixed_array[seq_index][offset] += 1
                        granularity_locks.release(seq_index)

                        if sam_nucleotide != ref_nucleotide:
                            mismatch += 1

                self.nucleotide_mismatch.value += mismatch
                sam_idx += cigar.Count
                current_position += cigar.Count

            elif cigar.Code == "D":

                offset = GenomeEvidence.nucleotide_offset["-"]
                for idx in range(cigar.Count):
                    # Shared memory is protected with a granularity lock
                    granularity_locks.acquire(seq_index)
                    evidence_fixed_array[current_position + idx][offset] += 1
                    granularity_locks.release(seq_index)

                self.deletion.value += 1
                current_position += cigar.Count

            elif cigar.Code == "I":

                insert_sequence = sam_record.Sequence[sam_idx: (sam_idx+cigar.Count)]
                insert_record = ReadSamFile.InsertTuple( Contigid=sam_record.Rname
                                                       , Sequenceidx=current_position
                                                       , Insertseq=insert_sequence)
                evidence_insert_list.append(insert_record)  # a list of inserts

                sam_idx += cigar.Count
                self.insertion.value += 1

            elif cigar.Code == "S":
                sam_idx += cigar.Count

            elif cigar.Code == "N":
                current_position += cigar.Count

        return

    def __decode_cigar(self, cigar):

        cigar_list = []
        for cigar_item in re.findall(self.cigar_regex, cigar):
            cigar_code = cigar_item[-1]
            cigar_count = int(cigar_item[:-1])
            cigar_list.append(ReadSamFile.CigarItem(Code=cigar_code, Count=cigar_count))

        return cigar_list

    def __queue_insert_records(self, insert_record_queue):

        for contig_id, contig_evidence in self.genome_evidence.items():

            evidence_insert_list = contig_evidence.Contiginsertlist

            for insertion in evidence_insert_list:
                insert_list = [insertion.Contigid, insertion.Sequenceidx, insertion.Insertseq]
                insert_record_queue.put(insert_list)
