#!/bin/bash

#	The regular expression that defines the fastq files
FILES=/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/gcpu/samples/GC032543_*/runs/run.160610.HiSeq2000.FCB.lane*/current/result/*.R1.fastq.gz 

for f in $FILES
do
  echo "Processing $f file..."
  # take action on each file. $f stores the current file name
  python trim_reads.py -i $f -o /uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Data/ATAC-Seq/160526.NextSeq.FCA/

done