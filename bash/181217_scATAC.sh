#!/bin/bash

#PBS -l nodes=1:ppn=2,mem=8gb,walltime=00:10:00

source activate scATAC

refGenome=/wehisan/general/user_managed/grpu_naik.s_2/Indexes/human38_withERCC
bedtoolsPath=/cm/shared/apps/bedtools/2.17.0/bin/
openChromatinBed=/uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/dbrown0/Bioinformatics/Resources/refGene_hg19_TSS.bed
tmpDir=/home/brown.d/Data/tmp

samtools fastq /home/brown.d/Data/ATAC/atac_a2.bam > /home/brown.d/Data/ATAC/atac_a2.fq