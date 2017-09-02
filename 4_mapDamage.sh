#!/bin/bash

## Script for running mapDamage to get deamination patterns and rescale BAM file quality
## Written by M. Nieves-Colon - 3/29/2017

## Usage on cluster: sbatch 4_mapDamage.sh

## Necessary modules
module load parallel/20140822
module load mapdamage/2.0.6

### Generate deamination plots and rescale bam files with mapDamage (assumes reference is indexed)
echo "***********************************************".
echo "**** Estimate damage and rescale bamfile ******"
mkdir 5_mapDamage
cd 5_mapDamage
cp ../4_MappingFiltering/*uniq.bam .
echo -e "Sample \t Rescaled?" > mapDamage.RescaledList.txt

# Using parallel because it cuts processing time by at least half
find *.bam | parallel -j4 -k 'mapDamage -i {} -r /home/mnievesc/ref-seqs/hg19-25chr-rCRS.fa -q --rescale'  # -quiet 

echo "***********************************************"
echo " "
echo " "
cd ..
