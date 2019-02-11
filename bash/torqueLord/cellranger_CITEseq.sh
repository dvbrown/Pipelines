#!/bin/bash

source activate tenX
cellrangerExec=/stornext/HPCScratch/home/brown.d/src/cellranger-3.0.2/cellranger
transcriptomeRef=/wehisan/general/user_managed/grpu_naik.s_2/TW/Ref_Genomes/refdata-cellranger-mm10-1.2.0
libraries=/stornext/HPCScratch/home/brown.d/src/CITEseq/librariesCITE.csv
featRef=/stornext/HPCScratch/home/brown.d/src/CITEseq/citeSeqReference_Dec2018.csv

$cellrangerExec count --id=dan_CITEseq \
                   --transcriptome=$transcriptomeRef \
                   --libraries=$libraries --feature-ref=$featRef \
                   --chemistry=SC3Pv2 \
                   --expect-cells=30000
