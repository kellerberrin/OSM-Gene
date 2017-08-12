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
# and setup a logger to receive analysis output.
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
            description="OSM_GENE Organism Genome Comparison Functionality")

        # --dir
        parser.add_argument("--dir", dest="workDirectory", default="./Work/",
                            help=("The work directory where log files and data files are found."
                                  " Use a Linux style directory specification with trailing forward"
                                  " slash '/' (default './Work/')."
                                  " Important - to run OSM_GENE this directory must exist, it will not be created."))

        # --gff
        parser.add_argument("--gff", dest="gffFile", default="ref.gff",
                            help=("The gff3 (not GFF2 or GTF) file that contains the genes, exons, etc for the"
                                  " chromosome(s)/contiguous region(s) of interest"))

        # --fasta
        parser.add_argument("--fasta", dest="fastaFile", default="ref.fasta",
                            help=("The fasta file that contains the reference nucleotide sequence for the chromsome(s) "
                                  "of interest"))

        # --parent
        parser.add_argument("--parent", dest="parentFile", default="noParent",
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

        # --contig
        parser.add_argument("--contig", dest="contigRegion", default="*",
                            help=("Define which contiguous DNA region (chromosome/mitochondria) to process."
                                  " Defaults to '*' for all contiguous regions."))

        # --log
        parser.add_argument("--log", dest="logFilename", default="OSM_GENE.log",
                            help=('Log file. Appends the log to any existing logs (default "OSM_QSAR.log").'
                                  'The log file always resides in the work directory.'))
        # --newlog
        parser.add_argument("--newlog", dest="newLogFilename", default="nonewlog", nargs='?',
                            help='Flush an existing log file (file name argument optional, default "OSM_GENE.log").'
                                 'The log file always resides in the work directory.')

        # "minimum_variant_counts" in Jeremy Horst's code
        parser.add_argument("--mutantcount", dest="minMutantCount", default=20,
                            help=("The minimum SAM/BAM coverage for the a single nucleotide analyzed in the Mutant"
                                  " genome."))

        # "minimum_variant_proportion" in Jeremy Horst's code
        parser.add_argument("--mutantprop", dest="minMutantProportion", default=0.7,
                            help=("The min proportion of a single nucleotide analyzed in the Mutant genome that is at "
                                  " variance from the parent (wild-type) genome"))

        # "maximum_wildtype_proportion" in Jeremy Horst's code
        parser.add_argument("--maxparentprop", dest="maxParentProportion", default=0.1,
                            help=("The max proportion of a single nucleotide analyzed in the Parent (wild-type) genome"
                                  "that is at variance from the reference (fasta) genome"))

        # "maximum_wildtype_variant_counts" in Jeremy Horst's code
        parser.add_argument("--maxparentcount", dest="maxParentCount", default=20,
                            help=("The maximum SAM/BAM count of a single nucleotide analyzed in the Parent (wild-type)"
                                  " genome that is at variance from the reference (fasta) genome"))

        # "minimum_wildtype_total_counts" in Jeremy Horst's code
        parser.add_argument("--minparentcount", dest="minParentCount", default=10,
                            help=("The minimum SAM/BAM count of a single nucleotide analyzed in the Parent (wild-type) "
                                  "genome that is to be compared to the Mutant genome"))

        parser.add_argument("--processes", dest="processCount", default=4,
                            help=("The number of CPU processes (not threads!) assigned to processing SAM genome data"))

        parser.add_argument("--queuesize", dest="queueSize", default=1000000,
                            help=("The maximum number of SAM records held in the inter-process record queue"))

        parser.add_argument("--lockgranularity", dest="lockGranularity", default=1000,
                            help=("The number of nucleotide positions per inter-process write lock (less is faster)"))

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
            ExecEnv.args.newLogFilename = os.path.join(ExecEnv.args.workDirectory,"OSM_GENE.log")
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

        description_text = ("OSM_Gene.py is a python script to find genetic differences in (optional)"
        " parent and mutant pairs by comparing genomic sequencing reads aligned to the same reference genome."
        " If the optional --parent flag is specified then the mutant genome is compared to the parent genome and"
        " the reference genome. The parent genome is also compared to the reference genome. If the --parent flag"
        " is not specified then the --mutant is compared to the reference FASTA genome only."
        " This program takes the genome FASTA file, the corresponding gene model in GFF3 (only) format, and the"
        " (optional) parent and mutant SAM files, respectively. Source FASTQ files should be filtered for quality"
        " and aligned to the FASTA reference genome with a tool such as BowTie2 or bwa, and output in SAM format.")

        usage_text = ("Usage:$python OSM_Gene.py --fasta <ref.fasta> --gff <ref.gff> --mutant <mutant.sam>")

        help_text = ("Help:$python OSM_Gene.py --help")

        return description_text + "\n\n" + usage_text + "\n\n" + help_text

