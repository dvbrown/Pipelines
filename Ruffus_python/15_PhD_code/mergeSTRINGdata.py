#!/usr/bin/python2.7

#try and merge data from String that has many source nodes and map expression data to it by gene name
import sys

mapFile = sys.argv[1]
expressionFile = sys.argv[2]
map = open(mapFile, 'r')
expression = open(expressionFile, 'r')

def makeBigList(file):
  newList = []
  currentGroup = []
  for line in file:
    line = line.split('\t')
    currentGroup.append(line[0:10])
    
    print currentGroup 

makeBigList(map)

