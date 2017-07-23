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

from __future__ import print_function,  division

from OSMGFFParser import ParseGFFSeq
from OSMGeneAnalysis import GeneAnalysis

# ==============================================================================#

class OSMGenomeComparison(object):

    def __init__(self, args, log):

        # Shallow copies of the runtime environment.
        self.log = log
        self.args = args

    def comparison(self):

            parsed_gff = ParseGFFSeq(self.args, self.log).parsed_structure
            parent_gene_analysis = GeneAnalysis(self.args, self.log, parsed_gff, self.args.parentFile)
            mutant_gene_analysis = GeneAnalysis(self.args, self.log, parsed_gff, self.args.mutantFile)

