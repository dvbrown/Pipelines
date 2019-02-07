#!/bin/bash

#PBS -l nodes=1:ppn=2,mem=8gb,walltime=00:10:00

source activate tenX

cellrangerExec=/home/brown.d/src/cellranger-3.0.2/cellranger
transcriptomeRef=/wehisan/general/user_managed/grpu_naik.s_2/Ref_Genomes/refdata-cellranger-mm10-1.2.0
libraries=/home/brown.d/src/CITEseq/librariesCITE.csv
featRef=/home/brown.d/src/CITEseq/citeSeqReference_Dec2018.csv

$cellrangerExec count --id=dan_CITEseq \
                   --transcriptome=$transcriptomeRef \
                   --libraries=$libraries --feature-ref=$featRef \
                   --expect-cells=1000
