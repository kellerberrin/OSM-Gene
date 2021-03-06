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

import time

# Import the runtime environment object.
from OSMExecEnv import ExecEnv, __version__

from OSMGFFParser import ParseGFFSeq
from OSMReadSam import ReadSamFile
from OSMSamEvidence import SamFileCPPLibrary
from OSMGeneAnalysis import GeneAnalysis


# ===================================================================================================
# The application object
# ====================================================================================================


class OSMGenomeComparison(object):

    def __init__(self, args, log):

        # Shallow copies of the runtime environment.
        self.log = log
        self.args = args

    def comparison(self):

        parsed_gff = ParseGFFSeq(self.log, self.args.fastaFile, self.args.gffFile).get_parsed_structure()

        mutant_evidence = self.__cpp_read_sam(parsed_gff, self.args.mutantFile)

        mutant_evidence.set_sam_filename("Mutant")

        mutant_all_snp = mutant_evidence.get_all_snp(self.args.minMutantCount, self.args.minMutantProportion)

        GeneAnalysis(self.log).print_snp_stats(mutant_all_snp, self.args.outputFile)

        mutant_cds_snp = mutant_all_snp.filter_cds_snp()

        if self.args.parentFile != "noParent":

            parent_evidence = self.__cpp_read_sam(parsed_gff, self.args.parentFile)

            parent_evidence.set_sam_filename("Parent")

            parent_all_snp = parent_evidence.get_all_snp(self.args.minParentCount, self.args.minParentProportion)

            GeneAnalysis(self.log).print_snp_stats(parent_all_snp, self.args.outputFile)

            parent_cds_snp = parent_all_snp.filter_cds_snp()

            novel_mutant_snp = mutant_cds_snp.difference(parent_cds_snp)

            GeneAnalysis(self.log).print_snp_stats(novel_mutant_snp, self.args.outputFile)

    def __read_sam(self, parsed_gff, sam_filename):

        sam_file_reader = ReadSamFile( self.args
                                     , self.log
                                     , parsed_gff
                                     , self.args.queueSize
                                     , self.args.lockGranularity
                                     , self.args.processCount
                                     , sam_filename)

        return sam_file_reader.get_evidence_object()

    def __cpp_read_sam(self, parsed_gff, sam_filename):

        sam_file_reader = SamFileCPPLibrary( self.args
                                            , self.log
                                            , parsed_gff
                                            , self.args.queueSize
                                            , self.args.lockGranularity
                                            , self.args.processCount
                                            , sam_filename)

        return sam_file_reader.get_evidence_object()




# ===================================================================================================
# The program mainline.
# ====================================================================================================

def main():

    try:

        ExecEnv()  # Setup the runtime environment.

        ExecEnv.log.info("############ OSM_GENE %s Start Comparison ###########", __version__)
        ExecEnv.log.info("Command Line: %s", ExecEnv.cmdLine)

        OSMGenomeComparison(ExecEnv.args, ExecEnv.log).comparison()  # Do the comparison.

        ExecEnv.log.info("Command Line: %s", ExecEnv.cmdLine)
        ExecEnv.log.info("Elapsed seconds CPU time %f (all processors, assumes no GPU).", time.clock())
        ExecEnv.log.info("############ OSM_GENE %s End Comparison ###########", __version__)

    except KeyboardInterrupt:

        ExecEnv.log.warning("\n\n")
        ExecEnv.log.warning("Control-C pressed. Program terminates. Open files may be in an unsafe state.")

    except IOError:

        ExecEnv.log.fatal("File error. Check directories, file names and permissions."
                          ' Check the default work directory "--dir" and --help".')
        ExecEnv.log.fatal("OSM_GENE exits.")

    except SystemExit:

        ExecEnv.log.info(ExecEnv.print_description())
        ExecEnv.log.info("OSM_GENE exits.")

    finally:

        clean_up = None  # Placeholder for any cleanup code.

if __name__ == "__main__":
    main()