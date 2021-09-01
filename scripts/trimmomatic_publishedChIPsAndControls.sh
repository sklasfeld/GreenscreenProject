#!/bin/bash

trimmomatic_install_dir="/usr/src"

# run Trimmomatic and remove TruSeq adapters
while read samp; do
    java -jar ${trimmomatic_install_dir}/Trimmomatic-0.39/trimmomatic-0.39.jar \
        SE -phred33 fastq/raw/${samp}.fastq.gz \
        fastq/trimmed/${samp}.trimmed.gz \
        ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33

    fastqc -o data/FASTQC/trimmed fastq/trimmed/${samp}.trimmed.fastq
done < meta/chip_controls_sampleIDs.list
