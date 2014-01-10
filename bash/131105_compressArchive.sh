#!/bin/bash

#compress the 2nd argument with file name 1st argument

tar -zcvf batch1_countExon.tgz GIC_011_out/GIC_011.countExons.txt GIC_020_out/GIC_020.countExons.txt \
GIC_034_out/GIC_034.countExons.txt GIC_035_out/GIC_035.countExons.txt \
GIC_039_out/GIC_039.countExons.txt GIC_041_out/GIC_041.countExons.txt
