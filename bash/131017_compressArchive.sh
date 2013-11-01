#!/bin/bash

#compress the 2nd argument with file name 1st argument

tar -zcvf batch1_HTSeqRev.tgz GIC_011/GIC_011_CGATGT.trim.bowtie.merged.sortS.sortP.dedup.sortName.revStrHTgem.txt  \
GIC_020/GIC_020_TGACCA.trim.bowtie.merged.sortS.sortP.dedup.sortName.revStrHTgem.txt GIC_034/GIC_034_ACAGTG.trim.bowtie.merged.sortS.sortP.dedup.sortName.revStrHTgem.txt \
GIC_035/GIC_035_GCCAAT.trim.bowtie.merged.sortS.sortP.dedup.sortName.revStrHTgem.txt GIC_039/GIC_039_CAGATC.trim.bowtie.merged.sortS.sortP.dedup.sortName.revStrHTgem.txt \
iGIC_041/GIC_041_CTTGTA.trim.bowtie.merged.sortS.sortP.dedup.sortName.revStrHTgem.txt 
