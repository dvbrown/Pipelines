#!/bin/bash

java -Xmx2g -jar /usr/local/picard/1.69/lib/AddOrReplaceReadGroups.jar \
INPUT=s_4_CD133n_B.bam.sorted.bam OUTPUT=./s_4_CD133n_B.bam.sorted.addReadGroupHeader.bam \
SORT_ORDER=unsorted RGID=s_4_CD133n_Blibrary \
RGLB=CD133n_B RGPL=illumina RGPU=FLOWCELL001 RGSM=CD133n;

java -Xmx2g -jar /usr/local/picard/1.69/lib/AddOrReplaceReadGroups.jar \
INPUT=s_4_CD133n_A.bam.sorted.bam OUTPUT=./s_4_CD133n_A.bam.sorted.addReadGroupHeader.bam \
SORT_ORDER=unsorted RGID=s_4_CD133n_Alibrary \
RGLB=CD133n_A RGPL=illumina RGPU=FLOWCELL001 RGSM=CD133n;

java -Xmx2g -jar /usr/local/picard/1.69/lib/AddOrReplaceReadGroups.jar \
INPUT=s_4_CD133n_C.bam.sorted.bam OUTPUT=./s_4_CD133n_C.bam.sorted.addReadGroupHeader.bam \
SORT_ORDER=unsorted RGID=s_4_CD133n_Clibrary \
RGLB=CD133n_C RGPL=illumina RGPU=FLOWCELL001 RGSM=CD133n;

java -Xmx2g -jar /usr/local/picard/1.69/lib/AddOrReplaceReadGroups.jar \
INPUT=s_5_CD133n_C.bam.sorted.bam OUTPUT=./s_5_CD133n_C.bam.sorted.addReadGroupHeader.bam \
SORT_ORDER=unsorted RGID=s_5_CD133n_Clibrary \
RGLB=CD133n_C RGPL=illumina RGPU=FLOWCELL001 RGSM=CD133n;

java -Xmx2g -jar /usr/local/picard/1.69/lib/AddOrReplaceReadGroups.jar \
INPUT=s_5_CD133n_A.bam.sorted.bam OUTPUT=./s_5_CD133n_A.bam.sorted.addReadGroupHeader.bam \
SORT_ORDER=unsorted RGID=s_5_CD133n_Alibrary \
RGLB=CD133n_A RGPL=illumina RGPU=FLOWCELL001 RGSM=sCD133n;

java -Xmx2g -jar /usr/local/picard/1.69/lib/AddOrReplaceReadGroups.jar \
INPUT=s_5_CD133n_B.bam.sorted.bam OUTPUT=./s_5_CD133n_B.bam.sorted.addReadGroupHeader.bam \
SORT_ORDER=unsorted RGID=s_5_CD133n_Blibrary \
RGLB=CD133n_B RGPL=illumina RGPU=FLOWCELL001 RGSM=CD133n;

java -Xmx2g -jar /usr/local/picard/1.69/lib/AddOrReplaceReadGroups.jar \
INPUT=s_6_CD133n_A.bam.sorted.bam OUTPUT=./s_6_CD133n_A.bam.sorted.addReadGroupHeader.bam \
SORT_ORDER=unsorted RGID=s_6_CD133n_Alibrary \
RGLB=CD133n_A RGPL=illumina RGPU=FLOWCELL001 RGSM=CD133n;

java -Xmx2g -jar /usr/local/picard/1.69/lib/AddOrReplaceReadGroups.jar \
INPUT=s_6_CD133n_B.bam.sorted.bam OUTPUT=./s_6_CD133_B.bam.sorted.addReadGroupHeader.bam \
SORT_ORDER=unsorted RGID=s_6_CD133n_Blibrary \
RGLB=CD133n_B RGPL=illumina RGPU=FLOWCELL001 RGSM=CD133n;

java -Xmx2g -jar /usr/local/picard/1.69/lib/AddOrReplaceReadGroups.jar \
INPUT=s_6_CD133n_C.bam.sorted.bam OUTPUT=./s_6_CD133n_C.bam.sorted.addReadGroupHeader.bam \
SORT_ORDER=unsorted RGID=s_6_CD133n_Clibrary \
RGLB=CD133n_C RGPL=illumina RGPU=FLOWCELL001 RGSM=CD133n;