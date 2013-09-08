#!/bin/bash

java -Xmx4g  -classpath /vlsci/VR0002/shared/Trimmomatic-0.22/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE -threads 1 -phred33 -trimlog $1'trimlog.txt' \
$1 $2 $1'pair.fq.gz' $1'unpair.gz' $2'pair.fq.gz' $2'unpair.gz' \
ILLUMINACLIP:/vlsci/VR0238/shared/pilotRawData/rawFastq/IlluminaAdapters.fa:2:40:15 LEADING:20 TRAILING:20 MINLEN:75
