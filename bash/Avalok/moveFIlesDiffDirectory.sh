#!/bin/bash
input="/path/to/txt/file"
while IFS= read -r var
do
  echo "$var"
  mv $var /uz/data/avalok/symbiosys/gcpi_r_kul_thierry_voet/projects/daniel_MEL006/bam_files/MEL06_T4_notSequenced
done < "$input"
