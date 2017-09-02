#!/bin/bash

## Script for removing adapters, trimming, and merging using SeqPrep.
## Written by T. Honap - June 26, 2015 and modified by M. Nieves Colon - Dec 1, 2015

## This script uses the program SeqPrep to remove residual adapter sequences from NGS reads.
## The program is run using default parameters. This can be changed depending on user requirements. 
## The -f and -r options specify the input R1 and R2 files respectively. These can be in .fastq.gz format. User can change the name as required. 
## The -1 and -2 options specify the output R1 and R2 files to be written after trimming. 
## The -s option specifies the final output file - the trimmed.merged.fastq.gz file. This file should be used for further analyses. 
## The -L option ensures that all trimmed and merged sequences having a length of less than 30 are discarded.
## The -A and -B options specify the forward and reverse adapters respectively. 

## All the sample folders must be in a main folder (Eg. HiSeqRun_June2015) and the script must also be placed in this main folder.
## Make sure the read files are named in the same format as mentioned here and below in the script (Eg. AD58_L001_R1_001.fastq.gz and AD58_L001_R2_001.fastq.gz). If necessary, modify the file names OR modify the script. 

## Usage on ASU Saguaro cluster: sbatch 1_SeqPrep.sh

## Necessary module
module load seqprep/jan2017

# To create list of file names
ls *R1_001.fastq.gz | cut -d "_" -f 1  >  list

cat list | while read line
 do
  echo "Processing sample $line"
  SeqPrep -f "${line}"_R1_001.fastq.gz -r "${line}"_R2_001.fastq.gz -o 11 -1 "${line}"_R1.trimmed.fastq.gz -2 "${line}"_R2.trimmed.fastq.gz  -3  "${line}"_R1.discarded.fastq.gz -4 "${line}"_R2.discarded.fastq.gz -L 30 -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -s "${line}".trimmed.merged.fastq.gz
 done


### References
# 1. http://www.linfo.org/wildcard.html
# 2. http://stackoverflow.com/questions/12174947/removing-a-part-of-filename-of-a-bunch-of-files
# 3. http://unix.ittoolbox.com/groups/technical-functional/unixadmin-l/unix-command-to-cut-some-part-of-the-filename-for-20-files-at-one-time-5150240
# 4. https://www.maketecheasier.com/rename-files-in-linux
