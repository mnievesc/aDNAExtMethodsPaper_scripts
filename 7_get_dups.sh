#!/bin/bash

## Script for getting number of PCR duplicates for clonality estimates with shotgun data
## Written by M. Nieves-Colon - 11/28/2017

## Usage: ./7_get_dups.sh

# Copy over mapped bam files
cp *trimmed.merged.mapped.bwa.bam . Recalculate_clonality/

# Sort and index BAM with samtools and the mark dups with Picard
for a in *bam;
 do 
  ID=$(basename $a) 
  NAME=$(echo $ID |cut -d "." -f1)   
  echo "Processing sample $NAME"    
  samtools sort $NAME.1Mds.trimmed.merged.mapped.bwa.bam $NAME.1Mds.trimmed.merged.mapped.bwa.sort
  samtools index $NAME.1Mds.trimmed.merged.mapped.bwa.sort.bam
  java -jar /usr/local/share/java/picard.jar MarkDuplicates INPUT=$NAME.1Mds.trimmed.merged.mapped.bwa.sort.bam OUTPUT=$NAME.1Mds.trimmed.merged.mapped.bwa.dedup.bam METRICS_FILE=$NAME.metrics.txt
done

# To get numbers
mkdir ../markdup-metrics_alldsmapped
mv *metrics.txt ../markdup-metrics_alldsmapped
