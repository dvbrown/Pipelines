# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 11:25:21 2016

@author: u0107775
"""
import pandas as pd
import os, re

table = '/Users/u0107775/Data/Mitochondria_Deletion/Fastq/Fastq_files/pindel/testPipeline/testConfig.filter.txt'
output = '/Users/u0107775/Data/Mitochondria_Deletion/Fastq/Fastq_files/pindel/testPipeline/testConfig.comma.txt'

#df = pd.read_csv(table, sep='\t|,')
df = pd.read_csv(table, sep='\t')

commaColumn = df.filter(like="AD")

# Split delimited values in a DataFrame column into two new columns
commaColumn = df['new_col1'], df['new_col2'] = zip(*df['sample2.AD'].apply(lambda x: x.split(',', 1)))

print commaColumn

# Make a new column 

#dfComma = df.join(s.apply(lambda x: Series(x.split(','))))

dfComma = df
dfComma.to_csv(output, sep='\t')
#os.system('head -n 25 /Users/u0107775/Data/Mitochondria_Deletion/Fastq/Fastq_files/pindel/testPipeline/testConfig.comma.csv')