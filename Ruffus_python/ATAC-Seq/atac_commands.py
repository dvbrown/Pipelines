import os, time
import pandas as pd

#################################    GLOBAL PARAMETERS    #####################################

# The reference genome is version 19 from Josien. The chr prefix is included
referenceGenomePath = '/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/jhaan0/humangenome/fasta/hg19Mt.fa'
binaryPath = '/cm/shared/apps/'

#################################    BEGIN COMMANDS    #####################################

def runJob(comm, taskName):
    '''An internal function used by the rest of the functions to spawn a process in the shell, 
    run a command based on the value of the string "comm" and capture the standard output.
    Throws an exception when failure occurs'''
    started = time.strftime('%X %x %Z')
    print '\n##############################################    RUNNNG TASK ' + taskName + ' at {0}'.format(started) +   '    ###############################################'
    print comm + '\n'
    #run the command
    os.system(comm)
    
def trimReads(inputFile, outputFile):
    'Take the raw sequencing reads and trim off the adpaters'
    comm = '''
    '''
    runJob(comm, 'TRIMMING READS')

def alignReads(inputFile, outputFile):
    'Align the fastq reads using bwa or bowtie or something'
    comm = '''
    '''
    runJob(comm, 'ALIGNING READS')