#!/bin/bash

#compress the 2nd argument with file name 1st argument

tar -zcvf batch1_fusionResults.tgz GIC_011_out/results* GIC_020_out/results* \
GIC_034_out/results* GIC_035_out/results* \
GIC_039_out/results* GIC_041_out/results*
