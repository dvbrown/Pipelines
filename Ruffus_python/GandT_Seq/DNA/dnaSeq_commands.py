import os, time, re

#################################    GLOBAL PARAMETERS    #####################################

refGenome = '/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Bioinformatics/Resources/hg19Mt'
# NEED a reference with the ERCC sequences!!!

binaryPath = '/cm/shared/apps/'
javaPath = '/cm/shared/apps/jdk/1.7.0/bin/java'
picardPath = '/cm/shared/apps/picard/current/'
gatkPath = '/cm/shared/apps/gatk/current/'
bedtoolsPath = '/cm/shared/apps/bedtools/2.17.0/bin/'
localBinaryPath = '/home/dbrown0/miniconda3/envs/py2bioinf/bin/'
tmpDir = ''

#################################    BEGIN COMMANDS    #####################################

def runJob(comm, taskName):
    '''An internal function used by the rest of the functions to spawn a process in the shell, 
    run a command based on the value of the string "comm" and capture the standard output.
    Throws an exception when failure occurs'''
    started = time.strftime('%X %x %Z')
    print('\n##############################################    RUNNNG TASK ' + taskName + ' at {0}'.format(started) +   '    ###############################################')
    print(comm + '\n')
    #run the command. Comment out the line below to print only the command and not run it.
    os.system(comm)
    
    
def trimReads(inputFile, outputFile):
    'Take the raw sequencing reads and trim off the adpaters'
    read2 = re.sub('.R1.fastq.gz', '.R2.fastq.gz', inputFile)
    outputFile2 = re.sub('.R1.', '.R2.', outputFile)
    #	Trim the Nextera adapter sequences
    comm = '''
    '''.format()
    runJob(comm, 'TRIMMING READS')
    

def alignReads(inputFile, outputFile):
    '''Align the fastq reads using bowtie'''
    #	Extract the read 2 filename
    read2 = re.sub('.R1.', '.R2.', inputFile)
    sampleName = inputFile[90:162]
    libraryName = inputFile[117:136]
    #   Build the read group information
    rgID = '{0} --rg SM:{1} --rg PL:ILLUMINA --rg LB:{2} '.format(inputFile, libraryName, sampleName)
    comm = '''
    '''.format()
    runJob(comm, 'ALIGNING READS')
    
    
def mergeBamPipeline(inputFileNames, outputFile):
    'As some samples are split over 4 lanes of sequencing on a NextSeq combine the aligned files'
    bam1 = inputFileNames[0]
    bam2 = inputFileNames[1]
    bam3 = inputFileNames[2]
    bam4 = inputFileNames[3]
    comm = '''java -Xmx5g -jar {0}MergeSamFiles.jar \
    INPUT={1} INPUT={2} INPUT={3} INPUT={4} \
    OUTPUT={5} SORT_ORDER=coordinate CREATE_INDEX=true \
    '''.format(picardPath, bam1, bam2, bam3, bam4, outputFile, tmpDir)
    runJob(comm, 'MERGING BAM FILES')
    
    
def sortSamtools(inputFile, outputFile):
    comm = "{0}samtools/current/samtools sort {1} > {2}".format(binaryPath, inputFile, outputFile)
    runJob(comm, "SORTING ALIGNMENTS")
    

def indexSamtools(inputFile):
    com = "{0}samtools/current/samtools index {1}".format(binaryPath, inputFile)
    # No output file therefore invoke os.system directly
    print('\n##############################################    RUNNNG TASK INDEX SAMTOOLS     ###############################################\n')
    os.system(com)
    
    
def estimateLibComplexity(inputFile, outputFile):
    'Estimate the library size using Picard tools. This should be done before duplicate removal'
    comm = '''java -Xmx5g -jar {0}EstimateLibraryComplexity.jar \
     I={1} \
     O={2}'''.format(picardPath, inputFile, outputFile)
    runJob(comm, 'ESTIMATING LIBRARY SIZE')
    
