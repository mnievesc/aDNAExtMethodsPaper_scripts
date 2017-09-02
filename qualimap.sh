#!/bin/bash

## Script for running qualimap to obtain coverage statistics
## Must run after rescaling quality with mapDamage. Run in same directory where rescaled bamfiles are stored.
## Written by M. Nieves-Colon - 3/29/2017

## Usage on cluster: sbatch qualimap.sh

### Load modules
module load qualimap/2.13 

### Generate mapping stats report with Qualimap    
echo -e "Total mtDNA reads \t Mean Read Depth \t StdDev Read Depth \t % ref-seq covered >1x \t % ref-seq covered >2x " > temp.qmstats
echo -e "Sample" > bamnames

ls *chrM.bam   > bamlist 

cat bamlist | while read line
 do
 sample=`echo $line | cut -d"." -f1`
 echo ${sample} >> bamnames 
 echo "Generating mapping stats for ${sample}"
 qualimap bamqc -bam ${line} -outdir ${sample}.QualimapStats -outformat pdf  # Do not output stdout to screen
 cd ${sample}.QualimapStats
 nreads=`grep "number of reads" genome_results.txt |cut -d "=" -f 2`
 meanrd=`grep "mean coverageData" genome_results.txt| cut -d "=" -f2`
 sdrd=`grep "std coverageData" genome_results.txt| cut -d "=" -f2`
 cov1=`grep "reference with a coverageData >= 1X" genome_results.txt | awk '{print $4}'`
 cov2=`grep "reference with a coverageData >= 2X" genome_results.txt | awk '{print $4}'`
 echo -e "$nreads \t $meanrd \t $sdrd \t $cov1 \t $cov2" >> ../temp.qmstats
 cd ../

done 

paste bamnames temp.qmstats  > QualimapStats.txt
rm temp.qmstats bamlist bamnames
rm *.bam
