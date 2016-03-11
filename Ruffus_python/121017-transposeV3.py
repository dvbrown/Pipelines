'''
Created on Aug 29, 2012

@author: d.brown6
'''
import sys
import csv

fileName = sys.argv[1]
fileOpen = open(fileName, 'r')
result = []
for i in fileOpen:
    item = i.strip('\n')
    item = i.split(',')
    result.append(i)
result =  filter(None, item)

#still have issues with the " and ' delimiters at the end of the line
stripped = []
for i in result:
    i = i.replace('\r', "', ")
    stripped.append(i)


output = csv.writer(open('output.csv', 'w'))
output.writerow(stripped)
print stripped