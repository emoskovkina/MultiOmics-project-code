#!/bin/bash

# Define the targets file
targetsFile="$@"

# Get all the sample names
sampleNames=$(cut -f1 $targetsFile)

for i in $sampleNames
do
  echo $i
  
  fastqFile=$(grep $i $targetsFile | cut -f2)
  fastqName="${fastqFile%.*}"
  echo $fastqName
  docker run --rm -v $PWD:/data ncbi/sra-tools:3.0.1 fasterq-dump --outdir /data $fastqName
  docker run -v /home:/home/ -w=$PWD --rm zlskidmore/kallisto:0.46.0 kallisto quant -i gencode.v39.transcriptome.idx -o $i --single -t 8 -l 250 -s 50 $fastqFile
done
