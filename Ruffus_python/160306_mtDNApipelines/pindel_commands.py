import time, os, subprocess, csv

# Global parameters
picardLoc = '/Users/u0107775/Bioinformatics/picard-tools-2.0.1/'
bioinformaticsDir = '/Users/u0107775/Bioinformatics/resources/'
#referenceGenome = '/Users/u0107775/Bioinformatics/resources/NC_012920.1.fa'
referenceGenome = '/Users/u0107775/Bioinformatics/resources/rCS.fa'

def runJob(comm, taskName, flagFile):
    '''An internal function used by the rest of the functions to spawn a process in the shell, capture the standard output 
    and generate a touch file. Runs the command in a shell and throws an exception when failure occurs'''
    started = time.strftime('%X %x %Z')
    print '\n##############################################    RUNNNG TASK ' + taskName + ' at {0}'.format(started) +   ' ###############################################'
    print comm + '\n'
    #run the command
    os.system(comm)
    finished = time.strftime('%X %x %Z')
    open(flagFile , 'w').write(finished)
    
def pindel():
    'Run the pindel tool for calling indels'
    pass

def pindel2vcf(inputFile, outFiles):
    'Convert the pindel output to vcf'
    output, flagFile = outFiles
    comm = '''/Users/u0107775/Bioinformatics/pindel/pindel2vcf -P {0} \
    -r {1} -R NC_012920.1 -d 31-OCT-2014 \
    -v {2} -G'''.format(inputFile, referenceGenome, output) # the name and date of the reference genome The G gives the genotypes in GATK format
    runJob(comm, "run pidel for deletions", flagFile)
    
def filterVcf():
    'filter the vcf file for high confidence variants'
    pass

def vcfToTable(inputFile, outFiles):
    'Convert a vcf file to a table for more human readability'
    output, flagFile = outFiles # Below I do want more fields than this - check the vcf file
    comm = '''java -jar /Users/u0107775/Bioinformatics/GenomeAnalysisTK.jar \
     -R {0} -T VariantsToTable \
     -V {1} \
     -F CHROM -F POS -F ID -F QUAL -F ALT -F END -F PE \
     -GF RC -GF RCL -GF RCR -GF DR -GF DV \
     -o {2}'''.format(referenceGenome, inputFile, output)
    runJob(comm, "convert a vcf file to a table", flagFile)
    
def calculateAlleleFreq():
    'Calculate the allele frequency based on reference and variant allele counts'