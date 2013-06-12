#!/bin/python

#hopefully a script that gets rid of commas

import sys

file = open(sys.argv[1], 'r')
for line in file:
	line.replace(",","")
	print file