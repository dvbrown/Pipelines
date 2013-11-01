#!/bin/bash
# Created by the VLSCI job script generator for SLURM on x86
# Tue Sep 03 2013 14:56:24 GMT+1000 (EST)

# Partition for the job:
#SBATCH -p main

# The name of the job:
#SBATCH --job-name="RNAseQC"

# Account to charge quota.
#SBATCH --account=VR0002

# Maximum number of CPU cores used by the job:
#SBATCH --ntasks=1

# The amount of memory in megabytes per process in the job:
#SBATCH --mem-per-cpu=16384

# Send yourself an email when the job:
# ends successfully
#SBATCH --mail-type=END

# Use this email address:
#SBATCH --mail-user=dvbrown@student.unimelb.edu.au

# The maximum running time of the job in days-hours:mins
#SBATCH --time=0-1:0

# Run the job from the directory where it was launched (default):
# The modules to load:
module load java/1.7.0_25
module load python-gcc/2.7.5
module load picard/1.96
module load rnaseqc/1.1.7
module load samtools-intel/0.1.19

# The job command(s):
./130918_ruffusPipeline.py -i GIC_011_CGATGT_L001_R1_001.fastq.trim.bowtie.bam -vv
