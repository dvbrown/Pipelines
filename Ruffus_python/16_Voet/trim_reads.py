#!/usr/bin/env python

import argparse, os, re

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Parser
parser = argparse.ArgumentParser(description="""With the help of a shell script descends into each directory containing 
raw sequencing reads separated by barcode and lane. Trims the reads using cutapadpt of the Nextera adapter sequence.""")
parser.add_argument('-i', '--inputFile', required=True, help='''The gzipped fastq file that will be trimmed''')
parser.add_argument('-o', '--outputDirectory', required=False, help='The directory to write all files to')
args = parser.parse_args()

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

# Assign the input specifed from the command line to a variable
inputFile = args.inputFile
outputDir = args.outputDirectory

def trimReads(input_file, output_dir):
    'Take the raw sequencing reads and trim off the adpaters. The output will all be in one directory'
    #   Get read 2 filename using string subsitiution
    read2 = re.sub('.R1.fastq.gz', '.R2.fastq.gz', input_file)
    #   Get the output filename by taking the end of the full path of the input filename
    output_file = input_file[-92:]
    # Build the path for the output directory by concatenating the output directory and output file
    output_file = output_dir + output_file
    output_file2 = re.sub('.R1.fastq.gz', '.R2.fastq.gz', output_file)
    
    comm = '''/home/dbrown0/miniconda3/envs/py2bioinf/bin/cutadapt -q 20,20 --minimum-length 35 \
    -a CTGTCTCTTATA -A CTGTCTCTTATA \
    -o {2} -p {3} \
    {0} {1} \
    '''.format(input_file, read2, output_file, output_file2)
    print('RUNNING TRIMMING \n' + comm)
    os.system(comm)

def main():   
    trimReads(inputFile, outputDir)
    
if __name__ == '__main__':
    main()