#!/bin/bash

# Join fasts together from the same Lane

cat DB_DCs_ADT_S1_L001_R1_001.fastq.gz DB_DCs_HTO_S2_L001_R1_001.fastq.gz > DB_DCs_CITE_S1_L001_R1_001.fastq.gz
cat DB_DCs_ADT_S1_L001_R2_001.fastq.gz DB_DCs_HTO_S2_L001_R2_001.fastq.gz > DB_DCs_CITE_S1_L001_R2_001.fastq.gz
cat DB_DCs_ADT_S1_L002_R1_001.fastq.gz DB_DCs_HTO_S2_L002_R1_001.fastq.gz > DB_DCs_CITE_S1_L002_R1_001.fastq.gz
cat DB_DCs_ADT_S1_L002_R2_001.fastq.gz DB_DCs_HTO_S2_L002_R2_001.fastq.gz > DB_DCs_CITE_S1_L002_R2_001.fastq.gz
cat DB_DCs_ADT_S1_L003_R1_001.fastq.gz DB_DCs_HTO_S2_L003_R1_001.fastq.gz > DB_DCs_CITE_S1_L003_R1_001.fastq.gz
cat DB_DCs_ADT_S1_L003_R2_001.fastq.gz DB_DCs_HTO_S2_L003_R2_001.fastq.gz > DB_DCs_CITE_S1_L003_R2_001.fastq.gz
cat DB_DCs_ADT_S1_L004_R1_001.fastq.gz DB_DCs_HTO_S2_L004_R1_001.fastq.gz > DB_DCs_CITE_S1_L004_R1_001.fastq.gz
cat DB_DCs_ADT_S1_L004_R2_001.fastq.gz DB_DCs_HTO_S2_L004_R2_001.fastq.gz > DB_DCs_CITE_S1_L004_R2_001.fastq.gz
