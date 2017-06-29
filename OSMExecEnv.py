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

import os
import sys
import argparse
import logging

__version__ = "0.1"

from OSMGenomeComp import OSMGenomeComparison

# ===================================================================================================
# A utility class to parse the program runtime arguments
# and setup a logger to receive classification output.
# ===================================================================================================

class ExecEnv(object):
    """Utility class to setup the runtime environment and logging"""

    # Static class variables.

    args = None
    log = None
    cmdLine = ""
    genome_compare = None

    def __init__(self):
        """Parse runtime arguments on object creation and maintain the runtime environment"""

        # Start a console logger to complain about any bad args (file logger defined below).

        file_log_format = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        console_log_format = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        ExecEnv.log = self.setup_logging(console_log_format)

        # Parse the runtime args
        parser = argparse.ArgumentParser(
            description="OSM_GENE Plasmodium Genome Functionality for Open Source Malaria")

        # --dir
        parser.add_argument("--dir", dest="workDirectory", default="./Work/",
                            help=('The work directory where log files and data files are found.'
                                  ' Use a Linux style directory specification with trailing forward'
                                  ' slash "/" (default "./Work/").'
                                  " Important - to run OSM_GENE this directory must exist, it will not be created."))

        # --gff
        parser.add_argument("--gff", dest="gffFile", default="model.gff",
                            help=("The gff3 file that contains the genes, exons, etc for the chromsome(s) of interest"))

        # --fasta
        parser.add_argument("--fasta", dest="fastaFile", default="model.fasta",
                            help=("The fasta file that contains the reference nucleotide sequence for the chromsome(s) "
                                  "of interest"))

        # --parent
        parser.add_argument("--parent", dest="parentFile", default="parent.sam",
                            help=("The SAM file was generated by an aligner (bowtie, bwa) that read FASTQ files"
                                  " generated by sequencing the parent organism"
                                  " and assigned the sequences to a position with respect to the"
                                  " reference genome specified by the --fasta flag"))

        # --mutant
        parser.add_argument("--mutant", dest="mutantFile", default="mutant.sam",
                            help=("The SAM file was generated by an aligner (bowtie, bwa) that read FASTQ files"
                                  " generated by sequencing the MUTANT organism"
                                  " and assigned the sequences to a position with respect to the"
                                  " reference genome specified by the --fasta flag"))

        # --log
        parser.add_argument("--log", dest="logFilename", default="OSM_QSAR.log",
                            help=('Log file. Appends the log to any existing logs (default "OSM_QSAR.log").'
                                  'The log file always resides in the work directory.'))
        # --newlog
        parser.add_argument("--newlog", dest="newLogFilename", default="nonewlog", nargs='?',
                            help='Flush an existing log file (file name argument optional, default "OSM_QSAR.log").'
                                 'The log file always resides in the work directory.')

      # --version
        parser.add_argument("--version", action="version", version=__version__)

        ExecEnv.args = parser.parse_args()

################# File Logging is here ########################################################

        ExecEnv.args.logFilename = os.path.join(ExecEnv.args.workDirectory,ExecEnv.args.logFilename)

        if ExecEnv.args.newLogFilename != "nonewlog" and ExecEnv.args.newLogFilename is not None:
            ExecEnv.args.newLogFilename = os.path.join(ExecEnv.args.workDirectory,ExecEnv.args.newLogFilename)
            log_append = False
            self.setup_file_logging(ExecEnv.args.newLogFilename, log_append, file_log_format)

        elif ExecEnv.args.newLogFilename != "nonewlog":  # No filename supplied (optional arg).
            ExecEnv.args.newLogFilename = os.path.join(ExecEnv.args.workDirectory,"OSM_QSAR.log")
            log_append = False
            self.setup_file_logging(ExecEnv.args.newLogFilename, log_append, file_log_format)

        else:
            log_append = True
            self.setup_file_logging(ExecEnv.args.logFilename, log_append, file_log_format)


        cmd_line = ""
        for argStr in sys.argv:
            cmd_line += argStr + " "

        ExecEnv.cmdLine = cmd_line


################# Other house keeping ########################################################

        # Check that the work directory exists and terminate if not.
        if not os.path.isdir(ExecEnv.args.workDirectory):
            ExecEnv.log.error('The OSM_GENE work directory: "%s" does not exist.', ExecEnv.args.workDirectory)
            ExecEnv.log.error("Create or Rename the work directory.")
            ExecEnv.log.error('Please examine the --dir" and "--help" flags.')
            sys.exit()

        ExecEnv.args.gffFile = os.path.join(ExecEnv.args.workDirectory, ExecEnv.args.gffFile)
        ExecEnv.args.fastaFile = os.path.join(ExecEnv.args.workDirectory, ExecEnv.args.fastaFile)
        ExecEnv.args.parentFile = os.path.join(ExecEnv.args.workDirectory, ExecEnv.args.parentFile)
        ExecEnv.args.mutantFile = os.path.join(ExecEnv.args.workDirectory, ExecEnv.args.mutantFile)

################# Create Comparison Instance ########################################################

        ExecEnv.genome_compare = OSMGenomeComparison(ExecEnv.args, ExecEnv.log)

################# File logging ########################################################

    def setup_logging(self, log_format):
        """Set up Python logging"""

        logger = logging.getLogger("OSMLogger")
        logger.setLevel(logging.INFO)  # Default output level.

        # Create a console log

        console_log = logging.StreamHandler()
        console_log.setLevel(logging.DEBUG)  # Output debug to screen
        console_log.setFormatter(log_format)

        logger.addHandler(console_log)

        return logger

    def setup_file_logging(self, log_filename, append, log_format):
        """Set up Python logging to log file"""

        # Create a file log.

        if append:
            file_log = logging.FileHandler(log_filename, mode='a')
        else:
            file_log = logging.FileHandler(log_filename, mode='w')

        file_log.setLevel(logging.INFO)  # Info level and above to file.
        file_log.setFormatter(log_format)

        ExecEnv.log.addHandler(file_log)
        if not append:
            ExecEnv.log.info("Flushed logfile: %s", log_filename)
        ExecEnv.log.info("Logging to file: %s", log_filename)

    @staticmethod
    def comparison():

        ExecEnv.genome_compare.comparison()

    @staticmethod
    def print_description():

        print("""
        MinorityReport.py is a python script meant to find genetic differences in parent-child diads or strain pairs by comparing genomic sequencing reads aligned to the same reference genome. 

        This program takes the genome FASTA file, the corresponding gene model in GFF3 format, and the parent and child SAM files, respectively. FASTQ files should be filtered for quality, aligned to a reference genome with a tool such as BowTie2, and output in SAM format. 

        The script MinorityReport-MASTER.py is meant to divide this task across processors on your computer, one per chromosome in the genome. This accelerates the process enormously, but requires one processor per chromosome.
        """, file=sys.stderr)
        print("\nUsage: MinorityReport.py <ref seq FASTA> <ref seq gff> <sam alignment parent> <sam alignment mutant>",
              file=sys.stderr)
        print(
            "\ne.g.:    ./MinorityReport.py species.fna species.gff species-bug1.bowtie2.sam species-bug2.bowtie2.sam",
            file=sys.stderr)
        print("\nNote: if using GTF or GFF2 format, be sure the file ends with '.gtf'", file=sys.stderr)
        print(
            "\nOptions: -vp  <minimum_variant_proportion>	    minimum fraction of reads at position supporting variant to accept (default=0.3).",
            file=sys.stderr)
        print(
            "         -wp  <maximum_variant_proportion>	    maximum fraction of reads at position supporting wildtype to accept (default=0.01).",
            file=sys.stderr)
        print(
            "         -vc  <minimum_variant_counts>	        minimum mutant-strain reads covering position to evaluate (default=30).",
            file=sys.stderr)
        print(
            "         -wc  <maximum_wildtype_variant_counts>	maximum parent-strain reads with variant to report (default=0).",
            file=sys.stderr)
        print(
            "         -wtc <minimum_wildtype_total_counts>	minimum total parent-strain reads covering position to report (default=0).",
            file=sys.stderr)
        print("", file=sys.stderr)
        print("         -cnv     analyze Copy Number Variants. Also available through CNV_caller.py", file=sys.stderr)
        print("         -keepsoft include read portions denoted as SoftClip for weak matches. default is to ignore.",
              file=sys.stderr)
        print("         -median  calculate median instead of mean CNV for tile range of positions (default: mean).",
              file=sys.stderr)
        print(
            "         -window_size       <sliding window length>  cnv: size of window to check for each position  (default=3000).",
            file=sys.stderr)
        print(
            "         -window_increment  <sliding window spacing> cnv: increment for progressing through chromosome(s) (default=3000).",
            file=sys.stderr)
        print(
            "         -report_frequency  <report window spacing>  cnv: increment for reporting - must be multiple of window_increment (default=3000).",
            file=sys.stderr)
        print("         -o_cnv   <output_file> specify output file for CNV analysis (default: print to stdout).",
              file=sys.stderr)
        print("", file=sys.stderr)
        print("         -pe  paired read matching alignments ONLY (default=False).", file=sys.stderr)
        print("         -gene_type  <GFF gene type>		descriptor to find correct entry lines in GFF file.",
              file=sys.stderr)
        print(
            "         -cds evaluate CDS for start and stop codons in each CDS, i.e. don't trust GFF file (default=False).",
            file=sys.stderr)
        #	print("         -all report all variants whether nonsynonymous or only nucleotide changes (default=False).",file=sys.stderr)
        print(
            "         -position_read_report <output_file> output a file with the number of counts that map to each position in the genome for parent & mutant (default=False).",
            file=sys.stderr)
        print("         -o   <output_file> specify output file (default: print to stdout).", file=sys.stderr)
        quit()


