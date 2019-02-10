#!/bin/bash

#PBS -l nodes=1:ppn=2,mem=8gb,walltime=00:10:00

source activate tenX

tagR1=/wehisan/home/allstaff/b/brown.d/SCORE/CITEseq/NN123/DB-DCs-ADT_S7_R1.fastq.gz
tagR2=/wehisan/home/allstaff/b/brown.d/SCORE/CITEseq/NN123/DB-DCs-ADT_S7_R1.fastq.gz
tagList=/home/brown.d/home/SCORE/CITEseq/citeSeq_Ab_tags.csv
outDir=/wehisan/home/allstaff/b/brown.d/SCORE/CITEseq/NN123/output

CITE-seq-Count -R1 $tagR1 -R2 $tagR2 -t $tagList \
-cbf 1 -cbl 16 -umif 17 -umil 26 -cells 30000 -o $outDir \
-n 1000
