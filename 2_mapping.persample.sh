#!/bin/bash

## Script for mapping and quslity filtering sequence reads using BWA and SAMtools.
## Written by M Nieves Colon - Dec 1, 2015

## Usage on cluster: sbatch mapping.persample.sh $NAME

## Reference sequence
ref=/home/mnievesc/ref-seqs/hg19-25chr-rCRS.fa 

## Necessary modules for Saguaro cluster
SAMTOOLS=/home/mnievesc/software/samtools-0.1.19/samtools
BWA=/home/mnievesc/software/bwa-0.7.5a/bwa 

## Working directory
cd /home/mnievesc/HiSeqrun2_Yale_March2017/


### 1. Map reads and filter with $SAMTOOLS. 
echo "****************************************"
echo "**** Mapping and filtering reads  ******"
mkdir 4_MappingFiltering   # change name to mapping and filtering
cd 4_MappingFiltering/redoing/

echo -e "Sample \t Mapped Reads \t  % Mapped \t Q30 Mapped Reads \t Removed Duplicates \t Duplicate Rate (Duplicates per Q30 Mapped reads) \t \
Mapped Unique Reads \t % endogenous (Mapped Unique/Total Reads) \t Mapped Reads average length \t \
Cluster Factor (Total mapped reads / Unique reads)" > $1.MappingFilteringStats.txt
cp ../../3_SeqPrep/$1.trimmed.merged.fastq.gz .


echo "Mapping data for $1"
  
$BWA aln -l 1000 -n 0.01 $ref $1.trimmed.merged.fastq.gz > $1.trimmed.merged.sai   # Map with seed disabled
$BWA samse $ref $1.trimmed.merged.sai $1.trimmed.merged.fastq.gz > $1.trimmed.merged.bwa.all.sam  # generate SAM alignment
$SAMTOOLS view -bSh $1.trimmed.merged.bwa.all.sam > $1.trimmed.merged.bwa.all.bam  # Generate bam file
$SAMTOOLS view -bh -F4 $1.trimmed.merged.bwa.all.bam > $1.trimmed.merged.mapped.bwa.bam  # Filter unmapped reads (F4 flag))
$SAMTOOLS view -bh -q 30 $1.trimmed.merged.mapped.bwa.bam > $1.trimmed.merged.mapped.q30.bwa.bam  # Filter Q30 quality 
$SAMTOOLS sort $1.trimmed.merged.mapped.q30.bwa.bam $1.trimmed.merged.mapped.q30.bwa.sort # Sort alignment by leftmost coordinate 
$SAMTOOLS rmdup -s $1.trimmed.merged.mapped.q30.bwa.sort.bam $1.trimmed.merged.mapped.q30.bwa.sort.rmdup.bam # Remove duplicates
# Remove reads with multiple mappings	 
$SAMTOOLS view -h $1.trimmed.merged.mapped.q30.bwa.sort.rmdup.bam | grep -v 'XT:A:R'| grep -v 'XA:Z' |grep -v 'XT:A:M' | awk '{if($0~/X1:i:0/||$0~/^@/  )print $0}' | $SAMTOOLS view -bS - > $1.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam 

# Print statistics per sample
mapped=`$SAMTOOLS view -c $1.trimmed.merged.mapped.bwa.bam`
q30=`$SAMTOOLS view -c $1.trimmed.merged.mapped.q30.bwa.sort.bam`
rmdup=`$SAMTOOLS view -c $1.trimmed.merged.mapped.q30.bwa.sort.rmdup.bam`
uniq=`$SAMTOOLS view -c $1.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam`
length=`$SAMTOOLS view $1.trimmed.merged.mapped.q30.bwa.sort.rmdup.uniq.bam | awk '{SUM+=length($10);DIV++}END{print SUM/DIV}'`  #fixed error in length estimation

echo -e "$1 \t $mapped \t `echo "calculate % mapped here"` \t $q30 \t $rmdup \t `echo "scale=4;($q30-$rmdup)/$rmdup" | bc` \t $uniq \t `echo "calculate endogenous here"` \t $length \t `echo "scale=4;$uniq/$q30" | bc` " >> $1.MappingFilteringStats.txt
echo -e "$1 \t $mapped \t `echo "calculate % mapped here"` \t $q30 \t $rmdup \t `echo "scale=4;($q30-$rmdup)/$rmdup" | bc` \t $uniq \t `echo "calculate endogenous here"` \t $length \t `echo "scale=4;$uniq/$q30" | bc` " >> MappingFilteringStats.txt
echo " "

# Remove merged file to declutter
rm $1.trimmed.merged.fastq.gz
