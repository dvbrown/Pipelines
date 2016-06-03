#!/bin/bash

#	The regular expression that defines the fastq files
FILES=/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Data/ATAC-Seq/160526.NextSeq.FCA/GC031715_AA*.160511_ATAC_Seq_Tn5.160526.NextSeq.FCA.lane1.gcap_dev.R1.fastq.gz

for f in $FILES
do
	#	Define each fastq file by lane
	
   lane2="${f/lane1/lane2}"
   lane3="${f/lane1/lane3}"
   lane4="${f/lane1/lane4}"
   echo "Processing $f"
   echo "Processing $lane2"
   echo "Processing $lane3"
   echo "Processing $lane4"
   # take action on each file. $f stores the current file name
   python2.7 atac_pipeline.py -i $f -i $lane2 -i $lane3 -i $lane4

done