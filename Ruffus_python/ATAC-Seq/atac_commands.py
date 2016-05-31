import os, time, re

#################################    GLOBAL PARAMETERS    #####################################

# The reference genome is version 19 from Josien. The chr prefix is included
refGenome = '/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/jhaan0/humangenome/fasta/hg19Mt.fa'
binaryPath = '/cm/shared/apps/'
picardPath = '/cm/shared/apps/'
tmpDir = ''

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
    comm = '''{0}cutadapt/1.4.1/
    '''.format(binaryPath)
    runJob(comm, 'TRIMMING READS')
    

def alignReads(inputFile, outputFile):
    '''Align the fastq reads using bwa or bowtie or something.  
    Paired-end reads were aligned to hg19 or mm10 using BOWTIE2 using the parameter –X2000 allowing fragments of up to 2 kb to align.
    Duplicates were removed and library size was estimated using PICARD tools'''
    read1 = ''
    read2 = ''
    rgID = ''
    comm = '''{0}bowtie2 --local -p 8 --rg-id {1} -x {2} -1 {3} -2 {4}
    | samtools view -bS -o '{5} -
    '''.format(binaryPath, rgID, refGenome, read1, read2, outputFile)
    runJob(comm, 'ALIGNING READS')
    
    
def mergeBams(inputFile, outputFile):
    'As some samples are split over multiple lanes of sequencing combine the aligned files'
    bam2 = re.sub('_L001_','_L002_', inputFile)
    bam3 = re.sub('_001.','_002.', inputFile)
    bam4 = re.sub('_001.','_002.', bam2)
    
    comm = '''java -Xmx10g -jar {0}MergeSamFiles.jar 
    INPUT= {1} + INPUT= {2} ' INPUT= {3} + INPUT= {4} OUTPUT= {5}
    CREATE_INDEX=true MAX_RECORDS_IN_RAM=750000 TMP_DIR={6}
    '''.format(picardPath, inputFile, bam2, bam3, bam4, outputFile, tmpDir)

    runJob(comm, 'MERGING BAM FILES')
    
    
def sortSamtools(inputFile, outputFile):
    comm = "samtools sort {0} > {1}".format(inputFile, outputFile)
    runJob(comm, "SORTING ALIGNMENTS")
    

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