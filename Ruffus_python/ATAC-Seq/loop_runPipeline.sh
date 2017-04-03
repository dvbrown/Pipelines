#!/bin/bash

#	The regular expression that defines the fastq files
FILES=/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/projects/daniel_atac/bam_files/GC038133*sorted.dedup.bam 

for f in $FILES
do
  echo "Processing $f file..."
  # take action on each file. $f stores the current file name
  sge batch -o $f.out -e $f.err python2.7 atac_pipeline.py -i $f

done