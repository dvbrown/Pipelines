import time, os, subprocess, csv

#################################    GLOBAL PARAMETERS    #####################################

picardPath = '/Users/u0107775/Bioinformatics/picard-tools-2.0.1/'
referenceGenomePath = '/Users/u0107775/Bioinformatics/resources/rCS.fa'
pindelPath = '/Users/u0107775/Bioinformatics/pindel/pindel'
gatkPath = '/Users/u0107775/Bioinformatics/GenomeAnalysisTK.jar'

#################################    BEGIN COMMANDS    #####################################

def runJob(comm, taskName, flagFile):
    '''An internal function used by the rest of the functions to spawn a process in the shell, 
    run a command based on the value of the string "comm" and capture the standard output.
    Also generates a touch file. Throws an exception when failure occurs'''
    started = time.strftime('%X %x %Z')
    print '\n##############################################    RUNNNG TASK ' + taskName + ' at {0}'.format(started) +   '    ###############################################'
    print comm + '\n'
    #run the command
    os.system(comm)
    finished = time.strftime('%X %x %Z')
    open(flagFile , 'w').write(finished)
    
def pindel(inputFile, outFiles, pindelConfig):
    'Run the pindel tool for calling indels. The pindel config file is specified at the commandline'
    output, flagFile = outFiles
    # Use string formatting to populate the the commands that will be run in the shell
    comm = '{0} -f {1} -o outputPindel -i {2}'.format(pindelPath, referenceGenomePath, pindelConfig)
    runJob(comm, 'pindel to detect variants', flagFile)

def pindel2vcf(inputFile, outFiles):
    'Convert the pindel output to vcf which contains all different types of variants'
    output, flagFile = outFiles
    # Provide the name and date of the reference genome The 'G' gives the genotypes in GATK format
    comm = '''/Users/u0107775/Bioinformatics/pindel/pindel2vcf -P {0} \
    -r {1} -R NC_012920.1 \
    -d 31-OCT-2014 \
    -v {2} -G'''.format(inputFile, referenceGenomePath, output)
    runJob(comm, 'pindel2vcf to parse raw pindel file', flagFile)
    
def subsetVcf(inputFile, outFiles):
    '''Filter the vcf file for high confidence variants'''
    output, flagFile = outFiles
    comm = ''' java -jar {0} \
   -R {1} \
   -T VariantFiltration \
   -V {2} -o {3} \
   --filterExpression "SVTYPE == DEL && SVLEN < -150" \
   --filterName "No_deletion"'''.format(gatkPath, referenceGenomePath, inputFile, output)
    runJob(comm, 'selectVariants from GATK to subset variants', flagFile)

def vcfToTable(inputFile, outFiles):
    'Convert a vcf file to a table for more human readability'
    output, flagFile = outFiles
    comm = '''java -jar {0} \
     -R {1} -T VariantsToTable \
     -V {2} \
     -F CHROM -F POS -F END ID -F REF -F ALT -F SVTYPE -F SVLEN \
     -GF GT -GF AD \
     -o {3}'''.format(gatkPath, referenceGenomePath, inputFile, output)
    runJob(comm, "convert a vcf file to a table", flagFile)
    
def calculateAlleleFreq():
    'Calculate the allele frequency based on reference and variant allele counts'
    # Write a python program here that uses pandas
    pass