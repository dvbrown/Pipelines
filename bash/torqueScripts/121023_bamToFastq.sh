#!/bin/bash
# Created by the VLSCI job script generator
# Tue Oct 23 2012 14:06:03 GMT+1100 (EST)

# Queue for the job:
#PBS -q batch terri

# Account to be charged
#PBS -A VR0002

# The name of the job:
#PBS -N BamToFq

# Maximum number of CPU cores used by the job:
#PBS -l procs=1

# The amount of memory in gigabytes per process in the job:
#PBS -l pvmem=2gb

# Send yourself an email when the job (e)nds successfully
#PBS -m e

# Use this email address:
#PBS -M dvbrown@student.unimelb.edu.au

# The maximum running time of the job in days:hours:mins:secs
#PBS -l walltime=0:1:0:0

# Run the job from the directory where it was launched:
cd $PBS_O_WORKDIR

# The modules to load:
module load picard/1.69
module load java/1.6.0_20
module load bpipe/0.9.5.3

# The job command(s):
bpipe run 121023_bamToFq.pipe ./*.bam