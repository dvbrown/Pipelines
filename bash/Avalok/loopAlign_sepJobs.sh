#!/bin/bash

#	The regular expression that defines the fastq files
#FILES=/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Data/GandT_Seq/Validation_160616/DNA/GC032370_CGAGGCTG-CTCCTTAC.160601_160601_GandT_Seq_test.160616.NextSeq.FCA.lane1.gcap_dev.R1.fastq.gz

FILES=/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Data/GandT_Seq/Validation_160616/DNA/GC032370_AAGAGGCA-AGGCTTAG*.R1.fastq.gz 

err=/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Data/GandT_Seq/Validation_160616/DNA/logs/err_
out=/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Data/GandT_Seq/Validation_160616/DNA/logs/out_

for f in $FILES
do
  #echo "Processing $f"
  # Extract only the filename instead of the whole path
  filename = ${$f:(-100)}
  
  # take action on each file. $f stores the current file name
  echo sge batch -e $err$filename.txt -o $out$filename.txt python dnaSeq_pipeline.py -i $f 

done