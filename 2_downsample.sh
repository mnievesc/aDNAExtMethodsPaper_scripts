#!/bin/bash

## Script to downsample reads from shotgun data
## Written by M. Nieves-Colon - 2/12/2017
## To downsample I chose the lowest number of merged reads within an individual pair.

# Following instructions from GitHUB README page for seqtk:
> Subsample 10000 read pairs from two large paired FASTQ files (remember to use the same random seed to keep pairing):
`seqtk sample -s100 read1.fq 10000 > sub1.fq`
`seqtk sample -s100 read2.fq 10000 > sub2.fq`

# I used the following commands on Saguaro cluster to subsample smallest number in a pair after SeqPrep merging and adapter trimming:

seqtk=/home/mnievesc/software/seqtk/seqtk

# Downsample lowest amount of reads in a pair. Use today's date as seed. Do each sample individually.

cd 4_SeqPrep

$seqtk sample -s 02122017 aGB7-D.trimmed.merged.fastq.gz 5192848 > ../10_Downsample_perind/aGB7-D.ds.trimmed.merged.fastq
$seqtk sample -s 02122017 aGB7-H.trimmed.merged.fastq.gz 5192848 > ../10_Downsample_perind/aGB7-H.ds.trimmed.merged.fastq

$seqtk sample -s 02122017 T-251-D.trimmed.merged.fastq.gz 1329660 > ../10_Downsample_perind/T-251-D.ds.trimmed.merged.fastq
$seqtk sample -s 02122017 T-251-H.trimmed.merged.fastq.gz 1329660 > ../10_Downsample_perind/T-251-H.ds.trimmed.merged.fastq

$seqtk sample -s 02122017 PI-388-D.trimmed.merged.fastq.gz 1007161 > ../10_Downsample_perind/PI-388-D.ds.trimmed.merged.fastq
$seqtk sample -s 02122017 PI-388-H.trimmed.merged.fastq.gz 1007161 > ../10_Downsample_perind/PI-388-H.ds.trimmed.merged.fastq

$seqtk sample -s 02122017 PI-67-D.trimmed.merged.fastq.gz 6876556 > ../10_Downsample_perind/PI-67-D.ds.trimmed.merged.fastq
$seqtk sample -s 02122017 PI-67-H.trimmed.merged.fastq.gz 6876556 > ../10_Downsample_perind/PI-67-H.ds.trimmed.merged.fastq

$seqtk sample -s 02122017 PC-E24-D.trimmed.merged.fastq.gz 2482634 > ../10_Downsample_perind/PC-E24-D.ds.trimmed.merged.fastq
$seqtk sample -s 02122017 PC-E24-H.trimmed.merged.fastq.gz 2482634 > ../10_Downsample_perind/PC-E24-H.ds.trimmed.merged.fastq

$seqtk sample -s 02122017 PC-117-D.trimmed.merged.fastq.gz 1050327 > ../10_Downsample_perind/PC-117-D.ds.trimmed.merged.fastq
$seqtk sample -s 02122017 PC-117-H.trimmed.merged.fastq.gz 1050327 > ../10_Downsample_perind/PC-117-H.ds.trimmed.merged.fastq
