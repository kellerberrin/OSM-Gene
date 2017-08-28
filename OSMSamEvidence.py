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



##############################################################################################
#
#   Base class object to hold the results of reading the SAM file
#
##############################################################################################


class SamFileEvidence(object):

    EvidenceFields = namedtuple("EvidenceFields", "Contigrecord Contigfixedarray Contiginsertarray")

    def __init__(self, log, sam_filename):

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


