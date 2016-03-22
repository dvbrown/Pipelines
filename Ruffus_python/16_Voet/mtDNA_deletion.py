import time, os, re, subprocess, csv
import pandas as pd

# Global parameters
picardLoc = '/Users/u0107775/Bioinformatics/picard-tools-2.0.1/'
bioinformaticsDir = '/Users/u0107775/Bioinformatics/resources/'
#referenceGenome = '/Users/u0107775/Bioinformatics/resources/NC_012920.1.fa'
referenceGenome = '/Users/u0107775/Bioinformatics/resources/rCS.fa'

def runJob(comm, taskName, flagFile):
    '''An internal function used by the rest of the functions to spawn a process in the shell, capture the standard output 
    and generate a touch file. Runs the command in a shell and throws an exception when failure occurs'''
    started = time.strftime('%X %x %Z')
    print '\n################################################### RUNNNG TASK ' + taskName + ' at {0}'.format(started) + ' ###############################################'
    print comm + '\n'
    #run the command
    os.system(comm)
    finished = time.strftime('%X %x %Z')
    open(flagFile , 'w').write(finished)
    
def convertUnalignedBam(inputFile, outFiles):
    'The mark illumina adapters need a bam file'
    read1 = inputFile
    read2 = inputFile.replace('_R1_', '_R2_')
    output, flagFile = outFiles
    comm = """java -Xmx2G -jar {0}picard.jar FastqToSam \
    FASTQ={1} FASTQ2={2} \
    OUTPUT={3} \
    READ_GROUP_NAME=test \
    SAMPLE_NAME=test \
    LIBRARY_NAME=Illumina \
    PLATFORM_UNIT=H0164ALXX140820.2 \
    PLATFORM=illumina \
    SEQUENCING_CENTER=Brussels""".format(picardLoc, read1, read2, output)
    runJob(comm, "convertunaligned", flagFile)
    
def markAdapters(unalignmedBam, outFiles):
    output, flagFile = outFiles
    comm= """java -Xmx2G -jar {0}picard.jar MarkIlluminaAdapters \
    I={1} O={2} \
    M=markilluminaadapters_metrics.txt \
    TMP_DIR=/Users/u0107775/Bioinformatics/temp
    """.format(picardLoc, unalignmedBam, output)
    runJob(comm, "markAdapters", flagFile)
    
def alignMtDNAGATK(unalignedBam, outFiles):
    'Align the sequencing reads against the mitochondrai genome. Add read group info later'
    output, flagFile = outFiles
    # Extract sample name for read group name
    sampleName = unalignedBam[:-4]
    readGroup = '@RG\tID:{}\tLB:Unknown\tPL:illumina\tCN:Brussels'.format(sampleName)
    comm = """java -Xmx2G -jar {0}picard.jar SamToFastq \
    I={1} \
    FASTQ=/dev/stdout CLIPPING_ATTRIBUTE='XT' \
    CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
    TMP_DIR=/Users/u0107775/Bioinformatics/temp | \
    bwa mem -M -t 3 -R {4} -p {3} /dev/stdin | \
    java -Xmx3G -jar {0}picard.jar MergeBamAlignment \
    ALIGNED_BAM=/dev/stdin UNMAPPED_BAM={1} \
    OUTPUT={2} R={3} \
    CREATE_INDEX=true ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
    TMP_DIR=/Users/u0107775/Bioinformatics/temp
    """.format(picardLoc, unalignedBam, output, referenceGenome, readGroup)
    runJob(comm, "cleanAndAlignBam", flagFile)
    
def alignMtDNA(inputFile, outFiles):
    'Align the sequencing reads against the mitochondrai genome. Add read group info later'
    read1 = inputFile
    read2 = inputFile.replace('_R1_', '_R2_')
    output, flagFile = outFiles
    sampleName = inputFile[:-4]  
    readGroup = '"@RG\tID:{0}\tSM:{0}\tPL:Illumina"'.format(sampleName)
    # Extract sample name for read group name
    comm = """bwa mem -M -t 3 -R {4} {2} {0} {1} > {3}
    """.format(read1, read2, referenceGenome, output, readGroup)
    runJob(comm, "cleanAndAlignBam", flagFile)
    
def delly(inputFile, outFiles):
    'Run the best deletion tools I can find which is delly'
    output, flagFile = outFiles
    # Set the path to the boost libriaires
    #os.system('export DYLD_LIBRARY_PATH=/Users/u0107775/Bioinformatics/delly/src/modular-boost/stage/lib')
    comm = "/Users/u0107775/Bioinformatics/delly/src/delly -t DEL -o {2} -g {1} {0}".format(inputFile, referenceGenome, output)
    runJob(comm, "run delly for deletions", flagFile)
    
def viewVcfFile(inputFile, outFiles):
    'Delly emits a binary file. Need to parse it to human readable'
    output, flagFile = outFiles
    comm = "bcftools view {} > {}".format(inputFile, output)
    runJob(comm, "translate bcf file", flagFile)
    
def indelList():
    'Produce a list of know indels'
    comm = """java -Xmx1g -jar {0}GenomeAnalysisTK.jar \
     -T RealignerTargetCreator \
     -R {0}rCRS_Mitochondira_fasta_noLines.fa \
     -o {0}output.intervals \
     -known {0}Mills_and_1000G_gold_standard.indels.hg38.vcf""".format(bioinformaticsDir)
    os.system(comm)

def localRealignment(bamFile, outFiles):
    'Perform local realignment'
    output, flagFile = outFiles
    comm = """ java -Xmx4g -Djava.io.tmpdir=/path/to/tmpdir \
    -jar {0}GenomeAnalysisTK.jar \
    -R  \
    -T IndelRealigner \
    -targetIntervals  \
    -o  \
    -known /path/to/indel_calls.vcf
    -LOD 0.4"""
    runJob(com, "extractSecondaryAlignments", flagFile)

def extractSecondaryAlignments(bamFile, outFiles):
    output, flagFile = outFiles
    os.system("samtools view -H {} > header.sam".format(bamFile))
    com = "samtools view {} | grep 'SA:' > {}".format(bamFile, output)
    runJob(com, "extractSecondaryAlignments", flagFile)

def addSamHeader(samFile, outFiles):
    output, flagFile = outFiles
    com = "cat header.sam {} | samtools view -bh > {}".format(samFile, output)
    runJob(com, "addSamHeader", flagFile)
    os.system("rm header.sam")

def sortSamtools(bamFile, outFiles):
    output, flagFile = outFiles
    com = "samtools sort {} > {}".format(bamFile, output)
    print com + "\n"
    runJob(com, "sortSamtools", flagFile)

def indexSamtools(bamFile, outFiles):
    output, flagFile = outFiles
    com = "samtools index {}".format(bamFile)
    # No touch file therefore invoke os.system directly
    os.system(com)
    
#def convertToBed(bamfile, outFiles):
#    '''Convert the bam file of chimeric reads to bed and retain the secondary 
#    alignment flag as the score field in the bed file'''
#    output, flagFile = outFiles
#    # Write a temporary bed file
#    os.system('bedtools bamtobed -i {} > temp.bed'.format(bamfile))
#    # Extract the SA tag from same bam file and paste it as a column to the temp bed file
#    com = "samtools view -h {0} | awk {{'print $17'}} | paste temp.bed - > {1}".format(bamfile, output)
#    runJob(com, 'BamToBed SA', flagFile)
#    # Clean up bed file
#    os.system('rm temp.bed')
    
def convertToBed(bamfile, outFiles):
    '''Sort the bam file by coordinate to have split reads next to each other
    The convert to bed'''
    output, flagFile = outFiles
    com = "samtools sort -n {0} | bedtools bamtobed -i > {1}".format(bamfile, output)
    runJob(com, 'BamToBed SA', flagFile)
    
def collapseSplitReads(bedFile, outFiles):
    'Read in a bedfile as a pandas dataframe then output split reads side by side in a table'
    output, flagFile = outFiles
    # Read in data and split the first occurance of a read by its duplicate entries (alternate mappings)
    df = pd.read_csv(bedFile, sep="\t", header=None)
    dfLeft = df[df.duplicated(subset=3, keep='first')]
    dfRight = df[~df.duplicated(subset=3, keep='first')]
    # Merge by read name
    #dfJoin = pd.merge(dfLeft, dfRight, how='inner', on=3)
    dfJoin = pd.merge(dfRight, dfLeft, how='inner', on=3)
    # Rename columns
    dfJoin.columns = ['chr_L', 'start_L', 'end_L', 'read_name', 'map_score_L', 'strand_L', 
    'chr_R', 'start_R', 'end_R', 'map_score_R', 'strand_R']
    # Set index of dataframe to read name, modifying the original object
    dfJoin.set_index('read_name', inplace=True)
    # Write to file
    dfJoin.to_csv(output, sep='\t')
    # Write the flagFile
    finished = time.strftime('%X %x %Z')
    open(flagFile , 'w').write(finished)
    
        
def filterSplitReads(bedFile, outFiles):
    'This script is work in progress - go up to '
    df = pd.read_csv("None-D99-1_S1_L001_bwa_RG_NA.sorted.SA.sort.split.bed", sep="\t", index_col = 0)

    # Filter reads that span the circular chromosome
    dfNo_origin = df[(df.start_L != 0) &  (df.start_R != 0)]
    
    # Filter split reads that align to separate strands
    dfStrand = dfNo_origin[dfNo_origin.strand_L == dfNo_origin.strand_R]
    
    dfStrand.to_csv('test.bed', sep='\t')
    # Next, only the deletions that are sequenced more than once for both sense and antisense strands 
    # and with identical breakpoints will be selected. The number of reads in both senses will be listed and summed.
    
def filterInsertSize():
    'This command filters the bam file for inserts larger than 497'
    comm = '''samtools view None-D99-1_S1_L001_bwa_RG_NA.sorted.bam | awk '{ if ($9 >= 497) print }' > testInsert3SD.sam'''