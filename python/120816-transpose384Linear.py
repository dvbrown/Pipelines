#!/usr/bin/python2.7
#A script to parse 384 well files and transpose to column vector

import sys

fileName = sys.argv[1]
file = open(fileName, 'U')

def TransposeTab(file):
#A function to transpose a tab-delimited text file and transpose to a column vector  
  result = []
  for line in file:
    line.strip('')
    line.split('\t')
    result.append(line +',')
  return result
  
transpose = TransposeTab(file)
transpose = str(transpose)
print transpose

outputFile = open('output.txt', 'w')
outputFile.write(transpose)

#outputFile.close()
file.close()