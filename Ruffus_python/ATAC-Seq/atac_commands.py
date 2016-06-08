import os, time, re

#################################    GLOBAL PARAMETERS    #####################################

#   The reference genome is version 19 from Josien. The chr prefix and mitochondrail DNA is included
refGenome = '/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Data/hg19Mt'
binaryPath = '/cm/shared/apps/'
javaPath = '/cm/shared/apps/jdk/1.7.0/bin/java'
picardPath = '/cm/shared/apps/picard/current/'
gatkPath = '/cm/shared/apps/gatk/current/'
bedtoolsPath = '/cm/shared/apps/bedtools/2.17.0/bin/'
localBinaryPath = '/home/dbrown0/.local/bin/'
nucleoAtacPath = '/home/dbrown0/local/NucleoATAC/bin/'
openChromatinBed = '/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Bioinformatics/Resources/GM.nucpos.bed'
tmpDir = '/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Data/ATAC-Seq/160526.NextSeq.FCA/tmp'

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
    read2 = re.sub('.R1.fastq.gz', '.R2.fastq.gz', inputFile)
    outputFile2 = re.sub('', '', outputFile)
    #	Trim the Nextera adapter sequences
    comm = '''{5}cutadapt -q 15,15 --minimum-length 35 \
    -a CTGTCTCTTATA -A CTGTCTCTTATA \
    -o {3} -p {4} {1} {2} \
    '''.format(binaryPath, inputFile, read2, outputFile, outputFile2, localBinaryPath)
    runJob(comm, 'TRIMMING READS')
    

def alignReads(inputFile, outputFile):
    '''Align the fastq reads using bowtie'''
    #	Extract the read 2 filename
    read2 = re.sub('.R1.fastq.gz', '.R2.fastq.gz', inputFile)
    sampleName = inputFile[90:162]
    libraryName = inputFile[117:136]
    #   Build the read group information
    rgID = '{0} --rg SM:{1} --rg PL:ILLUMINA --rg LB:{2} '.format(inputFile, libraryName, sampleName)
    comm = '''{0}bowtie2/current/bowtie2 --local -p 8 --rg-id {1} -x {2} -1 {3} -2 {4} \
    | {0}samtools/current/samtools view -bS -o {5} -S \
    '''.format(binaryPath, rgID, refGenome, inputFile, read2, outputFile)
    runJob(comm, 'ALIGNING READS')
    
    
def mergeBamPipeline(inputFileNames, outputFile):
    'As some samples are split over multiple lanes of sequencing combine the aligned files'
    bam1 = inputFileNames[0]
    bam2 = inputFileNames[1]
    bam3 = inputFileNames[2]
    bam4 = inputFileNames[3]
    comm = '''java -Xmx5g -jar {0}MergeSamFiles.jar \
    INPUT={1} INPUT={2} INPUT={3} INPUT={4} \
    OUTPUT={5} SORT_ORDER=coordinate \
    '''.format(picardPath, bam1, bam2, bam3, bam4, outputFile, tmpDir)
    runJob(comm, 'MERGING BAM FILES')
    
    
def sortSamtools(inputFile, outputFile):
    comm = "{0}samtools/current/samtools sort {1} > {2}".format(binaryPath, inputFile, outputFile)
    runJob(comm, "SORTING ALIGNMENTS")
    

def indexSamtools(inputFile):
    com = "{0}samtools/current/samtools index {1}".format(binaryPath, inputFile)
    # No output file therefore invoke os.system directly
    print '\n##############################################    RUNNNG TASK INDEX SAMTOOLS     ###############################################\n'
    os.system(com)
    
    
def removeDuplicates(inputFile, outputFile):
    'Remove duplicates using Picard'
    comm = '''java -Xmx5g -jar {0}MarkDuplicates.jar \
    INPUT={1} OUTPUT={2} METRICS_FILE={2}.txt \
    CREATE_INDEX=true \
    REMOVE_DUPLICATES=true \
    '''.format(picardPath, inputFile, outputFile)
    runJob(comm, 'REMOVING DUPLICATES')
    
    
def estimateLibComplexity(inputFile, outputFile):
    'Estimate the library size using Picard tools. This should be done before duplicate removal'
    comm = '''java -Xmx5g -jar {0}EstimateLibraryComplexity.jar \
     I={1} \
     O={2}'''.format(picardPath, inputFile, outputFile)
    runJob(comm, 'ESTIMATING LIBRARY SIZE')
    
def countAlignChr(inputFile, outputFile):
    'Obtain the number of reads mapping to each chromosome'
    comm = '''samtools idxstats {0} | cut -f 1,3 > {1}'''.format(inputFile, outputFile)
    runJob(comm, 'COUNT READS PER CHROMOSOME')


def removeMtDNAreads(inputFile, outputFile):
    'Remove mitochondrial reads from ATAC-Seq data'
    comm = '''samtools idxstats {0} | cut -f 1 | \
    grep -v chrM | xargs samtools view -b {0} > \
    {1}'''.format(inputFile, outputFile)
    runJob(comm, 'REMOVING MITOCHONDRIAL READS')
    
    
def nucleoatac(inputFile, outputFile):
    'Call nucleosomes using the nucleotac software by the Greenleaf lab: http://nucleoatac.readthedocs.io/en/latest/nucleoatac/'
    comm = '''{4}nucleoatac run --bed {0} \
    --bam {1} --fasta {2}.fa \
    --out {3} --write_all \
    '''.format(openChromatinBed, inputFile, refGenome, outputFile, nucleoAtacPath) 
    runJob(comm, 'RUNNING NUCLEOATAC')
    
# FRAGMENT ANALYSIS - from the single-cell ATAC-Seq paper
#As in our previous work3, we adjusted the plus strand aligning reads by +4 and the minus strand aligning reads by -5 bp 
#to represent the center of the transposon binding event. For calculating accessibility for each peak, 
#we counted the number of fragments (not reads) from each cell falling into each of the 50,000 peaks. 
#We therefore filtered libraries by requiring >15% of fragments falling in open chromatin (peak set defined above) and having a library size >10,000 
#as estimated by PICARD (Figure 1d).
