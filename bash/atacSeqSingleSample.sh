#!/bin/bash

# launch tophat2 aligner with bowtie2

refGenome="/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/TuxedoSuite_Ref_Files/Homo_sapiens\
/Ensembl/GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"
refTranscriptome="/vlsci/VR0002/shared/Reference_Files/Indexed_Ref_Genomes/TuxedoSuite_Ref_Files\
/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
picardPath=".jar"

samtools fastq $1> output.fastq

cutadapt -q 15,15 --minimum-length 35 \
  -a CTGTCTCTTATA -A CTGTCTCTTATA \
  -o output input

bowtie2 --local -p 8 --rg-id {1} -X 2000 -I 0 -x {2} -U {3} \
| {0}samtools/current/samtools view -bS -o {4} -S \
# binaryPath, rgID, refGenome, inputFile, outputFile

samtools sort {1} > {2}
#".format(binaryPath, inputFile, outputFile)

samtools index {1}
#.format(binaryPath, inputFile)

java -Xmx5g -jar {0}MarkDuplicates.jar \
  INPUT={1} OUTPUT={2} METRICS_FILE={2}.txt \
  CREATE_INDEX=true \
  REMOVE_DUPLICATES=true

java -Xmx5g -jar {0}EstimateLibraryComplexity.jar \
     I={1} \
     O={2}
#.format(picardPath, inputFile, outputFile)

samtools idxstats {0} | cut -f 1,3 > {1}
#.format(inputFile, outputFile)

samtools idxstats {0} | cut -f 1 | \
  grep -v chrM | xargs samtools view -b {0} > {1}

java -Xmx5g -jar {0}CollectInsertSizeMetrics.jar \
    I={1} O={2}.txt H={2}.pdf M=0.5
#format(picardPath, inputFile, outputFile)

python {0}krisDavie_makeHeatmap.py --sumsOnly -proc 8 {1} {2} 2000 {3}
#format(tssPipeup_exec, inputFile, openChromatinBed, outputFile)
