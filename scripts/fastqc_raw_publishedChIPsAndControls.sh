#!/bin/bash

# create output directory
mkdir -p data/FASTQC/raw

while read samp; do
    fastqc -o data/FASTQC/raw fastq/raw/${samp}.fastq.gz
done < meta/chip_controls_sampleIDs.list
