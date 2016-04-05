import os, time
import pandas as pd

#################################    GLOBAL PARAMETERS    #####################################

referenceGenomePath = '/Users/u0107775/Bioinformatics/resources/rCS.fa'
pindelPath = '/Users/u0107775/Bioinformatics/pindel/'
gatkPath = '/Users/u0107775/Bioinformatics/GenomeAnalysisTK.jar'

#################################    BEGIN COMMANDS    #####################################

def runJob(comm, taskName):
    '''An internal function used by the rest of the functions to spawn a process in the shell, 
    run a command based on the value of the string "comm" and capture the standard output.
    Also generates a touch file. Throws an exception when failure occurs'''
    started = time.strftime('%X %x %Z')
    print '\n##############################################    RUNNNG TASK ' + taskName + ' at {0}'.format(started) +   '    ###############################################'
    print comm + '\n'
    #run the command
    os.system(comm)
    
def pindel(pindelConfig, outputFile):
    'Run the pindel tool for calling indels. The pindel config file is specified at the commandline'
    # Use string formatting to populate the the commands that will be run in the shell
    comm = '{0}pindel -f {1} -o {2} -i {3}'.format(pindelPath, referenceGenomePath, outputFile, pindelConfig)
    runJob(comm, 'PINDEL')
    # Write and empty file to signify that this step successfully completed. This is ncessary because pindel does not generate an output file
    f = open(outputFile, 'w')
    f.close()

def pindel2vcf(inputFile, outputFile):
    'Convert the pindel output to vcf which contains all different types of variants'
    # Provide the name and date of the reference genome The 'G' gives the genotypes in GATK format
    comm = '''{0}pindel2vcf \
    -P {1} \
    -r {2} \
    -R NC_012920.1 \
    -d 31-OCT-2014 \
    -v {3} -G'''.format(pindelPath, inputFile, referenceGenomePath, outputFile)
    runJob(comm, 'PINDEL TO VCF')
    
def subsetVcf(inputFile, outputFile):
    '''Filter the vcf file for high confidence variants.
    Select only deletions and deletions greater than 150bp.'''
    comm = ''' java -jar {0} \
   -R {1} \
   -T VariantFiltration \
   -V {2} -o {3} \
   --filterExpression "SVTYPE != 'DEL' || SVLEN > -150" \
   --filterName "Not_Deleted \
   invertFilterExpression"'''.format(gatkPath, referenceGenomePath, inputFile, outputFile)
    runJob(comm, 'SELECT VARIANTS')

def vcfToTable(inputFile, outputFile):
    'Convert a vcf file to a table for more human readability'
    comm = '''java -jar {0} \
     -R {1} -T VariantsToTable \
     -V {2} \
     -F CHROM -F POS -F END -F REF -F ALT -F SVTYPE -F SVLEN \
     -GF AD \
     -o {3}'''.format(gatkPath, referenceGenomePath, inputFile, outputFile)
    runJob(comm, "CONVERT vcf TO TABLE")
    
def calculateAlleleFreq(inputFile, outputFile):
    '''Calculate the allele frequency based on reference and variant allele counts'''
    # Read in the vcf to table file using both the tab and comma delimiters to split the allele counts into 2 columns
    df = pd.read_csv(inputFile, sep='\t|,')
    # Pandas only counts the number of columns that are tab delimited. Therefore write immediately to file and reload table later
    df.to_csv(outputFile, sep='\t')
    
    # Extract the general column names and then extract the sample specific column names
    colHeader = list(df.columns.values)[0:7]
    sampleNames = list(df.columns.values)[7:]
    
    # Initialize a new column header and for each sample append reference and deletion names
    sampleCounts = []
    for name in sampleNames:
        newName = name[:-3]
        refCounts = newName + '.REF'
        delCounts = newName + '.DEL'
        sampleCounts.append(refCounts)
        sampleCounts.append(delCounts)
    
    # Cncatenate the general and sample-specific column names
    newHeader = colHeader + sampleCounts
    
    # Read in the freshly written file that should have the correct number of columns.
    df2 = pd.read_csv(outputFile, sep='\t', skiprows=1, header=None, index_col=False)
    # Replace with the new column headers and write to file
    df2.columns = newHeader
    
    # Calculate the allele frequency
    for column in df2:
        # Extract columns containing reference read counts
        if '.REF' in column:
            # Extract sample names
            sampleName = column[:-4]
            # Generate column headers for indexing based on sample name
            refReads = '{}.REF'.format(sampleName)
            delReads = '{}.DEL'.format(sampleName)
            # Convert the deletions to floating point number, otherwise integer division is performed
            numerator = df2[delReads].apply(float)
            # Write new column label
            columnLabel = sampleName + '.FREQ'
            # Perform allele frequency calculation
            df2[columnLabel] = (numerator / (numerator + df2[refReads])) * 100    
    
    df2.to_csv(outputFile, sep='\t')