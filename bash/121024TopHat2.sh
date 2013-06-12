#!/bin/bash

# launch tophat2 aligner with bowtie2

refGenome="/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/TuxedoSuite_Ref_Files/Homo_sapiens\
/Ensembl/GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"

refTranscriptome="/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/TuxedoSuite_Ref_Files\
/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"

tophat -p 8 -o $1'_out' --fusion-search --mate-inner-dist 70 \
-G $refTranscriptome --transcriptome-index=transcriptome_data/known \
$refGenome $2 $3