#!/bin/bash

cellrangerExec=/stornext/HPCScratch/home/brown.d/src/cellranger-3.0.2/cellranger
transcriptomeRef=/wehisan/general/user_managed/grpu_naik.s_2/TW/Ref_Genomes/refdata-cellranger-mm10-1.2.0
libraries=/stornext/HPCScratch/home/brown.d/Data/CITEseq/

$cellrangerExec count --id=danCITE_RNA_only \
  --transcriptome=$transcriptomeRef \
  --fastqs=$libraries \
  --sample=DB_DCs_RNA \
  --expect-cells=30000
