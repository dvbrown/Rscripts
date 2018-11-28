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
import sys, os, re
import atac_commands 

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
                        help="""Write some help""")     
#    parser.add_option("-o", "--output_directory", dest="output_directory",
#                        action="append",
#                        default = list(),
#                        metavar="FILE",
#                        type="string",
#                        help="""The directory where the output should go to.""")
                        
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

#   End preamble, begin pipeline 
#   The code below is from Josien sge batch
#   sge batch -e errorfiles/qc1 -o outputfiles/qc1 -N QC1 sh shfiles/qc1.sh\
#   output /uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Data/ATAC-Seq

#################################    PIPELINE CODE GOES HERE    #####################################

# Assign the input specifed from the command line to a variable
inputFile = options.input_file
#outputDir = options.output_directory

# def trimReads(input_file, output_dir):
#     'Take the raw sequencing reads and trim off the adpaters. The output will all be in one directory'
#     #   Get read 2 filename using string subsitiution
#     read2 = re.sub('.R1.fastq.gz', '.R2.fastq.gz', input_file)
#     #   Get the output filename by taking the end of the full path of the input filename
#     output_file = input_file[-92:]
#     # Build the path for the output directory by concatenating the output directory and output file
#     output_file = output_dir + output_file
#     output_file2 = re.sub('.R1.fastq.gz', '.R2.fastq.gz', output_file)
#     comm = '''/home/dbrown0/.local/bin/cutadapt -q 20,20 --minimum-length 35 \
#     -a CTGTCTCTTATA -A CTGTCTCTTATA \
#     -o {2} -p {3} \
#     {0} {1} \
#     '''.format(input_file, read2, output_file, output_file2)
#     print(comm)
#     os.system(comm)
    
#   trimReads(inputFile[0], outputDir[0])

#@transform(inputFile, suffix('.fastq.gz'), '.bam')
#def runAlignment(inputFile, outputFile):
#    atac_commands.alignReads(inputFile, outputFile)
#    
@transform(inputFile, suffix('Aligned.sortedByCoord.out.bam'), '.merge.bam')
def runBamMerge(inputFile, outputFile):
   atac_commands.mergeBams(inputFile, outputFile)
    
#   Align the fastqs from each lane in a single script
# @transform(inputFile, suffix('lane1.gcap_dev.R1.fastq.gz'), 'lane1.gcap_dev.R1.bam')
# def runAligLane1(inputFile, outputFile):
#     atac_commands.alignReads(inputFile, outputFile)
# 
# # Make this a follows decorator
# @follows(runAligLane1) 
# @transform(inputFile, suffix('lane2.gcap_dev.R1.fastq.gz'), 'lane2.gcap_dev.R1.bam')
# def runAligLane2(inputFile, outputFile):
#     atac_commands.alignReads(inputFile, outputFile)
# 
# @follows(runAligLane2) 
# @transform(inputFile, suffix('lane3.gcap_dev.R1.fastq.gz'), 'lane3.gcap_dev.R1.bam')
# def runAligLane3(inputFile, outputFile):
#     atac_commands.alignReads(inputFile, outputFile)
# 
# @follows(runAligLane3)
# @transform(inputFile, suffix('lane4.gcap_dev.R1.fastq.gz'), 'lane4.gcap_dev.R1.bam')
# def runAligLane4(inputFile, outputFile):
#     atac_commands.alignReads(inputFile, outputFile)
# 
# #   Merge all bams from different lanes together into one file
# #   Generate output file name. Make this using grep next time
# mergeName = inputFile[0]
# mergeName = mergeName[91:116]
# 
# @follows(runAligLane4)  
# @merge([runAligLane1, runAligLane2, runAligLane3, runAligLane4], '{0}.merge.bam'.format(mergeName))
# def runBamMergePipeline(inputFileNames, outputFile):
#     atac_commands.mergeBamPipeline(inputFileNames, outputFile)
    
#-------------------------    POST ALIGNMENT    -----------------------------

#@follows(runIndexing)
#@transform(inputFile, suffix('.bam'), '.lib_metrics.txt')
#def runEstimateLibraryStats(inputFile, outputFile):
#    atac_commands.estimateLibComplexity(inputFile, outputFile)
#    
#@follows(runEstimateLibraryStats)
#@transform(inputFile, suffix('.bam'), '.deDup.bam')
#def runDuplicateRemoval(inputFile, outputFile):
#    atac_commands.removeDuplicates(inputFile, outputFile)
#    
#@transform(inputFile, suffix('.bam'), '.chrs.txt')
#def runChrAlignStats(inputFile, outputFile):
#    atac_commands.countAlignChr(inputFile, outputFile)
#    
#@follows(runChrAlignStats)
#@transform(inputFile, suffix('.bam'), '.noMt.bam')
#def runMTremoval(inputFile, outputFile):
#    atac_commands.removeMtDNAreads(inputFile, outputFile)
#    
#@transform(runMTremoval, suffix('.bam'), '')
##	Samtools index does not generate an output using the standard output
#def runIndexing(inputFile, outputFile):
#    atac_commands.indexSamtools(inputFile)

#@transform(inputFile, suffix('.bam'), '')
#def runKrisATAC(inputFile, outputFile):
#    atac_commands.kDavieATAC(inputFile, outputFile)
    
# @transform(inputFile, suffix('.bam'), '')
# def runPicardInsert(inputFile, outputFile):
#     atac_commands.collectInsertSize(inputFile, outputFile)
    
#################################    END PIPELINE    #####################################

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
