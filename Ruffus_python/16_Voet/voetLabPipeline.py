#!/usr/bin/env python
"""

    ruffus_template.py  [--input_file]
                        [--log_file PATH]
                        [--verbose]
                        [--target_tasks]
                        [--jobs]
                        [--just_print]
                        [--flowchart]
                        [--key_legend_in_graph]
                        [--forced_tasks]

"""
import sys, os
# import tasks
# import postAlign
import mtDNA_deletion
import indelTools

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    from optparse import OptionParser
    import StringIO

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %progs [options]")

    #
    #   general options: verbosity / logging
    #
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="count", default=1,
                      help="Print more verbose messages for each additional verbose level.")
    parser.add_option("-L", "--log_file", dest="log_file",
                      metavar="FILE",
                      type="string",
                      help="Name and path of log file")

    #
    #   pipeline options
    #
    parser.add_option("-i", "--input_file", dest="input_file",
                        action="append",
                        default = list(),
                        metavar="FILE",
                        type="string",
                        help="""The file(s) to use as input. If there are multiple
                        files use the -i argument multiple times""")
    parser.add_option("-t", "--target_tasks", dest="target_tasks",
                        action="append",
                        default = list(),
                        metavar="JOBNAME",
                        type="string",
                        help="Target task(s) of pipeline.")
    parser.add_option("-j", "--jobs", dest="jobs",
                        default=1,
                        metavar="N",
                        type="int",
                        help="Allow N jobs (commands) to run simultaneously.")
    parser.add_option("-n", "--just_print", dest="just_print",
                        action="store_true", default=False,
                        help="Don't actually run any commands; just print the pipeline.")
    parser.add_option("--flowchart", dest="flowchart",
                        metavar="FILE",
                        type="string",
                        help="Don't actually run any commands; just print the pipeline "
                             "as a flowchart.")
    #
    #   Less common pipeline options
    #
    parser.add_option("--key_legend_in_graph", dest="key_legend_in_graph",
                        action="store_true", default=False,
                        help="Print out legend and key for dependency graph.")
    parser.add_option("--forced_tasks", dest="forced_tasks",
                        action="append",
                        default = list(),
                        metavar="JOBNAME",
                        type="string",
                        help="Pipeline task(s) which will be included even if they are up to date.")

    # get help string
    f =StringIO.StringIO()
    parser.print_help(f)
    helpstr = f.getvalue()
    (options, remaining_args) = parser.parse_args()


    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    #
    #   Add names of mandatory options,
    #       strings corresponding to the "dest" parameter
    #       in the options defined above
    #
    mandatory_options = [ ]

    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    def check_mandatory_options (options, mandatory_options, helpstr):
        """
        Check if specified mandatory options have b een defined
        """
        missing_options = []
        for o in mandatory_options:
            if not getattr(options, o):
                missing_options.append("--" + o)

        if not len(missing_options):
            return

        raise Exception("Missing mandatory parameter%s: %s.\n\n%s\n\n" %
                        ("s" if len(missing_options) > 1 else "",
                         ", ".join(missing_options),
                         helpstr))
    check_mandatory_options (options, mandatory_options, helpstr)


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   imports


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from ruffus import *
import subprocess

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions

def run_cmd(cmd_str):
    """
    Throw exception if run command fails
    """
    process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, shell = True)
    stdout_str, stderr_str = process.communicate()
    if process.returncode != 0:
        raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                            (cmd_str, stdout_str, stderr_str, process.returncode))

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Logger

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

import logging
logger = logging.getLogger(options.target_tasks)
#
# We are interesting in all messages
#
if options.verbose:
    logger.setLevel(logging.DEBUG)
    stderrhandler = logging.StreamHandler(sys.stderr)
    stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
    stderrhandler.setLevel(logging.DEBUG)
    logger.addHandler(stderrhandler)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


#   Pipeline


#################################    PIPELINE CODE GOES HERE    #####################################

inputFile = options.input_file

#@transform(inputFile, suffix(".fastq.gz"), ['.bam', 'align.txt'])
#def convertFq(inputFile, output):
#    'Align with BBmap'
#    bbmap.alignMtDNA(inputFile, output)

#@transform(inputFile, suffix(".fastq.gz"), ['.U.bam', 'unaligned.txt'])
#def convertFq(inputFile, output):
#    'Convert the fastq files to unaligned bam'
#    mtDNA_deletion.convertUnalignedBam(inputFile, output)
#
@transform(inputFile, suffix(".U.bam"), ['.txt', 'adapters.txt'])
def markIlumminaAdapters(inputFile, output):
    'Convert the fastq files to unaligned bam'
    mtDNA_deletion.markAdapters(inputFile[0], output)

@transform(inputFile, suffix(".fastq.gz"), ['.m.bam', 'mt.txt'])
def alignMtDNA(inputFile, output):
    'Align the cleaned files to the mitochondira DNA genome'
    mtDNA_deletion.alignMtDNA(inputFile, output)

@transform(alignMtDNA, suffix(".bam"), ['.s.bam', ".sort.txt"])
def sortBam(inputFile, output):
    'Sort reads with samtools'
    mtDNA_deletion.sortSamtools(inputFile[0], output)

@transform(sortBam, suffix(".bam"), ["", "index.txt"])
def index(inputFile, output):
    'Sort reads with samtools'
    mtDNA_deletion.indexSamtools(inputFile[0], output)
#
#@transform(sortBam, suffix(".bam"), ['.bcf', 'del.txt'])
#def runDelly(inputFile, output):
#    'Align the cleaned files to the mitochondira DNA genome'
#    mtDNA_deletion.delly(inputFile[0], output)
#
#@transform(runDelly, suffix(".bcf"), ['.vcf', 'view.txt'])
#def convertBcf(inputFile, output):
#    'Convert the bcf file to be human readable'
#    mtDNA_deletion.viewVcfFile(inputFile[0], output)

#@transform(inputFile, suffix(".bam"), ['.SA.sam', ".extractSA.txt"])
#def getSecAlignment(inputFile, outFiles):
#    'Extract reads with secondary alignment flag in Sam'
#    mtDNA_deletion.extractSecondaryAlignments(inputFile, outFiles)
#
#@transform(getSecAlignment, suffix(".sam"), ['.bam', ".catHeader.txt"])
#def catHeader(inputFile, outFiles):
#    'Extract reads with secondary alignment flag in Sam'
#    # Because the output of the previous step is an array of 2 take only the first element (ie not the touchfile)
#    mtDNA_deletion.addSamHeader(inputFile[0], outFiles)

#@transform(inputFile, suffix(".bam"), ['.bed', ".bam2Bed.txt"])
#def convertBed(inputFile, outFiles):
#    'Convert bam to Bed'
#    mtDNA_deletion.convertToBed(inputFile, outFiles)

#@transform(inputFile, suffix(".bed"), ['.split.bed', 'splitRead.txt'])
#def castSplitReads(inputFile, output):
#    'Extract the secondary split read from bed file and append as columns next to the first split read entry'
#    mtDNA_deletion.collapseSplitReads(inputFile, output)

#@transform(inputFile, suffix(".vcf"), ['.csv', "success.txt"])
#def runVcfToTable(inputFile, output):
#    'extract fields from the VCF to a table format that is more convenient to work with in downstream analyses.'
#    indelTools.vcfToTable(inputFile, output)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Print list of tasks

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if options.just_print:
    pipeline_printout(sys.stdout, options.target_tasks, verbose=options.verbose)


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Print flowchart

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
elif options.flowchart:
    # use file extension for output format
    output_format = os.path.splitext(options.flowchart)[1][1:]
    pipeline_printout_graph (open(options.flowchart, "w"),
                             output_format,
                             options.target_tasks,
                             no_key_legend = True)
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Run Pipeline

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
else:
    pipeline_run(options.target_tasks,  multiprocess = options.jobs,
                        logger = logger, verbose=options.verbose)