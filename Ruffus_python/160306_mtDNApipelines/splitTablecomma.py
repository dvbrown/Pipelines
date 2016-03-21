# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 11:25:21 2016

@author: u0107775
"""
import pandas as pd

table = '/Users/u0107775/Data/Mitochondria_Deletion/Fastq/Fastq_files/pindel/testPipeline/testConfig.filter.txt'
output = '/Users/u0107775/Data/Mitochondria_Deletion/Fastq/Fastq_files/pindel/testPipeline/testConfig.comma.txt'

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

df2.to_csv(output, sep='\t')