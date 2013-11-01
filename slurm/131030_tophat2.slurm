#!/bin/bash
# Created by the VLSCI job script generator for SLURM on x86
# Fri Sep 13 2013 16:28:47 GMT+1000 (EST)

# Partition for the job:
#SBATCH -p main

# SMP job: want exclusive use of 1 node
#SBATCH --nodes=1
#SBATCH --exclusive

# The name of the job:
#SBATCH --job-name="tophat"

# Account to charge quota.
#SBATCH --account=VR0002

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=32768

# Send yourself an email when the job:
# ends successfully
#SBATCH --mail-type=END

# Use this email address:
#SBATCH --mail-user=dvbrown@student.unimelb.edu.au

# The maximum running time of the job in days-hours:mins
#SBATCH --time=0-1:0

# Run the job from the directory where it was launched (default):
# The modules to load:
module load bowtie2-intel/2.1.0
module load samtools-intel/0.1.19
module load python-gcc/2.7.5
module load tophat/2.0.4

# The job command(s):
./alignPipeline.py -i GIC_011_CGATGT_L001_R1_001.fastq.trim.fastq -vvv
