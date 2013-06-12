#!/bin/bash

casavaPath= #path to CASAVA
exptPath= #path to experiment folder where HiSeq data is stored

$casavaPath/CASAVA/bin/configureBclToFastq.pl --input-dir $exptPath/Data/Intensities/BaseCalls --output-dir Unaligned \
--mismatches 1 --sample-sheet ./121016sampleSheetDanBrown.csv

cd Unaligned
nohup make -j 4