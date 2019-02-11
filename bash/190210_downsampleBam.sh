#!/bin/bash

inputBam=/
outputBam=/

frac=$( samtools idxstats $inputBam | cut -f3 | awk 'BEGIN {total=0} {total += $1} \
END {frac=15000000/total; if (frac > 1) {print 1} else {print frac}}' )

samtools view -bs $frac $inputBam > $outputBam
