#!/bin/bash

# run Trimmomatic and remove TruSeq adapters
while read samp; do
    java -jar ~/Downloads/Trimmomatic-0.36/trimmomatic-0.36.jar \
        SE -phred33 fastq/raw/${samp}.fastq.gz \
        fastq/trimmed/${samp}.trimmed.gz \
        ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33

    fastqc -o FASTQC/trimmed fastq/trimmed/${samp}.trimmed.fastq
done < meta/chip_controls_sampleIDs.list
