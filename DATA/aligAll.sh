#!/bin/bash

DIR_FILE=PairedEnd
bwa=bwa
ref=../../human_g1k_v37.fasta

for line in $(ls PairedEnd/*read1.fastq|sed 's/.read1.fastq//g'|awk -F"/" '{print $2}')
do

temp=$(echo $line)

echo "--------------------------------------------"
echo "Aligning $temp"

$bwa mem  $ref ${DIR_FILE}/$temp.read1.fastq ${DIR_FILE}/$temp.read2.fastq > ${DIR_FILE}/$temp.sam

samtools view -S -b ${DIR_FILE}/$temp.sam > ${DIR_FILE}/$temp.bam

samtools sort ${DIR_FILE}/$temp.bam ${DIR_FILE}/$temp.sorted

samtools index ${DIR_FILE}/$temp.sorted.bam

rm ${DIR_FILE}/$temp.bam

echo "-----------------------------------------------------------------"
echo "===============Complete $line"
echo '-----------------------------------------------------------------'

done 
