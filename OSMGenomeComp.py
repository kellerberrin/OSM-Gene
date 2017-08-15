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
from OSMGeneEvidence import ReadSamFile, GenomeEvidence
from OSMGeneAnalysis import GeneAnalysis
from OSM_Filter import GenomeSNPFilter

# ==============================================================================#

class OSMGenomeComparison(object):

    def __init__(self, args, log):

        # Shallow copies of the runtime environment.
        self.log = log
        self.args = args

    def comparison(self):

        parsed_gff = ParseGFFSeq(self.log, self.args.fastaFile, self.args.gffFile).get_parsed_structure()

        mutant_evidence = ReadSamFile( self.log
                                      , parsed_gff
                                      , self.args.queueSize
                                      , self.args.lockGranularity
                                      , self.args.processCount
                                      , self.args.mutantFile).get_evidence_object()

        snp_filter = GenomeSNPFilter(self.log, mutant_evidence, self.args.minMutantCount, self.args.minMutantProportion)
        snp_filter.filter_gene_snp()
        snp_filter.filter_cds_snp()

#        GeneAnalysis(self.log, mutant_evidence).print_contig_stats()
