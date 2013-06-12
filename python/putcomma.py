#/bin/python

#put in commas in whitespace

import sys

file = open(sys.argv[1], 'r')
for line in file:
	line.replace(' ', '\n')
	print line
	
