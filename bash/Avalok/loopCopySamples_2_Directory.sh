#!/bin/bash

# The output directory change this when you want a different directory
outputDir=/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Data/GandT_Seq/Validation_160616/RNA/

#	Read the file that contains barcodes
while read f; do

  echo copying /uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/gcpu/samples/GC032370_$f/runs/run.160616.NextSeq.FCA.lane*/current/result/*.R1.fastq.gz $outputDir

done <dna_samples.txt
# Change the sample file when you want to read in different barcodes. The sample file is just a list of barcodes one per line