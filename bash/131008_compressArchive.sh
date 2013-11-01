#!/bin/bash

#compress the 2nd argument with file name 1st argument

tar -zcvf 131008_trimMergeBam.tgz $2 *.trim.bowtie.merged.sortS.sortP.dedup.bam*
