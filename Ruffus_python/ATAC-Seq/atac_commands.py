import os, time
import pandas as pd

#################################    GLOBAL PARAMETERS    #####################################

# The reference genome is version 19 from Josien. The chr prefix is included
referenceGenomePath = '/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/jhaan0/humangenome/fasta/hg19Mt.fa'
binaryPath = '/cm/shared/apps/'

# Get this fastqc file /uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/gcpu/samples/GC031537_AAGAGGCA-GTAAGGAG/runs/run.160506.HiSeq2000.FCB.lane2/current/fastqc/current/result/GC031537_AAGAGGCA-GTAAGGAG.160428_Tn5_HCC38.160506.HiSeq2000.FCB.lane2.gcap_16_04.R1_fastqc

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
    '''Align the fastq reads using bwa or bowtie or something.  
    Paired-end reads were aligned to hg19 or mm10 using BOWTIE2 using the parameter –X2000 allowing fragments of up to 2 kb to align.
    Duplicates were removed and library size was estimated using PICARD tools'''
    comm = '''
    '''
    
    headParams = 'bowtie2 --local -p 8 --rg-id ' + rgID
    midParams = ' -x ' + refGenome + ' -1 ' + read1 + ' -2 ' + read2
    tailParams = ' | samtools view -bS -o ' + output + ' -'
    comm = headParams + midParams + tailParams    
    
    runJob(comm, 'ALIGNING READS')
    
    
def mergeBams(bamFile, outFiles):
    'As some samples are split over multiple lanes of sequencing combine the aligned files'
    comm = '''
    '''
    
    output, flagFile = outFiles
    sample2 = re.sub('_L001_','_L002_', bamFile)
    sample3 = re.sub('_001.','_002.', bamFile)
    sample4 = re.sub('_001.','_002.', sample2)
    #------------------------------build shell command--------------------------------------
    headParams = 'java -Xmx10g -jar /usr/local/picard/1.96/lib/MergeSamFiles.jar'
    midParams = ' INPUT=' + bamFile + ' INPUT=' + sample2 + ' INPUT=' + sample3 + ' INPUT=' + sample4
    tailParams = ' CREATE_INDEX=true MAX_RECORDS_IN_RAM=750000 TMP_DIR=/vlsci/VR0238/shared/tmp'
    comm = headParams + midParams + ' OUTPUT=' + output + tailParams
    #---------------------------------------------------------------------------------------
    
    runJob(comm, 'MERGING BAM FILES')
    
    
def sortSamtools(inputFile, outputFile):
    com = "samtools sort {} > {}".format(bamFile, output)
    print com + "\n"
    runJob(com, "SORTING ALIGNMENTS")
    

def indexSamtools(inputFile):
    com = "samtools index {}".format(inputFile)
    # No output file therefore invoke os.system directly
    os.system(com)
    
    
def removeDuplicates(inputFile, outputFile):
    'Remove duplicates using Picard'
    comm = '''
    '''
    runJob(comm, 'REMOVING DUPLICATES')
    
    
def estimateLibSize(inputFile, outputFile):
    'Estimate the library size using Picard tools.'
    comm = '''
    '''
    runJob(comm, 'ESTIMATING LIBRARY SIZE')
    
# FRAGMENT ANALYSIS - from the single-cell ATAC-Seq paper
#As in our previous work3, we adjusted the plus strand aligning reads by +4 and the minus strand aligning reads by -5 bp 
#to represent the center of the transposon binding event. For calculating accessibility for each peak, 
#we counted the number of fragments (not reads) from each cell falling into each of the 50,000 peaks. 
#We therefore filtered libraries by requiring >15% of fragments falling in open chromatin (peak set defined above) and having a library size >10,000 
#as estimated by PICARD (Figure 1d).