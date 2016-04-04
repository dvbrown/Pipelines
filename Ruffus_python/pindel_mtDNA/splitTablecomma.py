# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 11:25:21 2016

@author: u0107775
"""
import pandas as pd

table = '/Users/u0107775/Data/Mitochondria_Deletion/Fastq/Fastq_files/pindel/testPipeline/testConfig.filter.txt'
output = '/Users/u0107775/Data/Mitochondria_Deletion/Fastq/Fastq_files/pindel/testPipeline/testConfig.calc.txt'

# Read in the vcf to table file using both the tab and comma delimiters to split the allele counts
df = pd.read_csv(table, sep='\t|,')
# Pandas only counts the number of columns that are tab delimited. Therefore write immediately to file and reload table
df.to_csv(output, sep='\t')

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

# Read in the fresh file that should have the correct number of columns.
df2 = pd.read_csv(output, sep='\t', skiprows=1, header=None, index_col=False)
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
        denominator = (df2[refReads].apply(float)) + numerator
        # Write new column label
        columnLabel = sampleName + '_freq'
        # Perform allele frequency calculation
        df2[columnLabel] = (numerator / (numerator + df2[refReads])) * 100
        
print df2    
df2.to_csv(output, sep='\t')