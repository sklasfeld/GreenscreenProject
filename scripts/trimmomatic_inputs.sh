#!/bin/bash

trimmomatic_install_dir="/usr/src"

# run Trimmomatic without adapter removal
no_adapter_input_list=("B" "C" "D" "E"
        "F" "G" "H" "I" "K"
        "L" "N" "P" "Q" "T")
for x in "${no_adapter_input_list[@]}"; do
    java -jar ${trimmomatic_install_dir}/Trimmomatic-0.39/trimmomatic-0.39.jar \
        SE -phred33 fastq/raw/input${x}.fastq \
        fastq/trimmed/input${x}.trimmed.fastq \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33

    fastqc -o data/FASTQC/trimmed fastq/trimmed/input${x}.trimmed.fastq
done

# run Trimmomatic and remove TruSeq adapters
adapter_input_list=("A" "J" "K" "M" "O"
        "R" "S")
for x in "${adapter_input_list[@]}"; do
    java -jar ${trimmomatic_install_dir}/Trimmomatic-0.39/trimmomatic-0.39.jar \
        SE -phred33 fastq/raw/input${x}.fastq \
        fastq/trimmed/input${x}.trimmed.fastq \
        ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33

    fastqc -o data/FASTQC/trimmed fastq/trimmed/input${x}.trimmed.fastq
done
