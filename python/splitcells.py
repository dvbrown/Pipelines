#A script to split lines by commas on a spreadsheet.

import sys


argument = open(sys.argv[1], 'r')
file = argument.readlines()

for line in file:
	splited = line.split(",")
	stripped = line.strip('"')
	print stripped
	

	