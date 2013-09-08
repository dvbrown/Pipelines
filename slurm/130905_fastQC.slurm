#!/bin/bash
# Created by the VLSCI job script generator for SLURM on x86
# Tue Sep 03 2013 14:56:24 GMT+1000 (EST)

# Partition for the job:
#SBATCH -p main

# The name of the job:
#SBATCH --job-name="130903_fastQc"

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
module load fastqc/0.10.1
module load java/1.7.0_25

# The job command(s):
./130903_fastqc.sh \ 
Sample_GIC_011/GIC_011_CGATGT_L001_R2_001.fastq.gz Sample_GIC_020/GIC_020_TGACCA_L001_R2_001.fastq.gz \
Sample_GIC_034/GIC_034_ACAGTG_L001_R2_001.fastq.gz Sample_GIC_035/GIC_035_GCCAAT_L001_R2_001.fastq.gz \
Sample_GIC_039/GIC_039_CAGATC_L001_R2_001.fastq.gz Sample_GIC_041/GIC_041_CTTGTA_L001_R2_001.fastq.gz
