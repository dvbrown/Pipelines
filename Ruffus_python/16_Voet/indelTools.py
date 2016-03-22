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
    
def delly(inputFile, outFiles):
    'Run the best deletion tools I can find which is delly'
    output, flagFile = outFiles
    # Set the path to the boost libriaires
    #os.system('export DYLD_LIBRARY_PATH=/Users/u0107775/Bioinformatics/delly/src/modular-boost/stage/lib')
    comm = "/Users/u0107775/Bioinformatics/delly/src/delly -t DEL -o {2} -g {1} {0}".format(inputFile, referenceGenome, output)
    runJob(comm, "run delly for deletions", flagFile)
    
def pindel(inputFile, outFiles):
    'Run the pindel SV tool'
    output, flagFile = outFiles
    comm = '''/Users/u0107775/Bioinformatics/pindel/pindel2vcf -P {0} \
    -r {1} -R NC_012920.1 -d 31-OCT-2014 \
    -v {2} -G'''.format(inputFile, referenceGenome, output) # the name and date of the reference genome The G gives the genotypes in GATK format
    runJob(comm, "run pidel for deletions", flagFile)
    
def platypus(inputFile, outFiles):
    'run the platypus tool to find deletions'
    output, flagFile = outFiles
    comm = '''python /Users/u0107775/Bioinformatics/Platypus_0.8.1/Platypus.py callVariants \
    --bamFiles={0} --refFile={1} --output={2} reigons=chrM'''.format(inputFile, referenceGenome, output)
    runJob(comm, "run Platypus", flagFile)
    
def viewVcfFile(inputFile, outFiles):
    'Delly emits a binary file. Need to parse it to human readable'
    output, flagFile = outFiles
    comm = "bcftools view {} > {}".format(inputFile, output)
    runJob(comm, "translate bcf file", flagFile)
    
def vcfToTable(inputFile, outFiles):
    'Convert a vcf file to a table for more human readability'
    output, flagFile = outFiles
    comm = '''java -jar /Users/u0107775/Bioinformatics/GenomeAnalysisTK.jar \
     -R {0} -T VariantsToTable \
     -V {1} \
     -F CHROM -F POS -F ID -F QUAL -F ALT -F END -F PE \
     -GF RC -GF RCL -GF RCR -GF DR -GF DV \
     -o {2}'''.format(referenceGenome, inputFile, output)
    runJob(comm, "convert a vcf file to a table", flagFile)
    
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
    
def filterInsertSize():
    'This command filters the bam file for inserts larger than 497'
    comm = '''samtools view None-D99-1_S1_L001_bwa_RG_NA.sorted.bam | awk '{ if ($9 >= 497) print }' > testInsert3SD.sam'''