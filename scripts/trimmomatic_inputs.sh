#!/bin/bash

# hard-code directories
fastq_raw_dir="fastq/raw"
fastq_trim_dir="fastq/trimmed"
fastqc_trim_dir="data/FASTQC/trimmed"

# create output directories
mkdir -p ${fastq_trim_dir}
mkdir -p ${fastqc_trim_dir}

trimmomatic_install_dir="/usr/src"

# run Trimmomatic without adapter removal
no_adapter_input_list=("B" "C" "D" "E"
        "F" "G" "H" "I"
        "L" "N" "P" "Q" "T")
for x in "${no_adapter_input_list[@]}"; do
    raw_fastq="${fastq_raw_dir}/input${x}.fastq.gz"
    trim_fastq="${fastq_trim_dir}/input${x}.trimmed.fastq.gz"
    java -jar ${trimmomatic_install_dir}/Trimmomatic-0.39/trimmomatic-0.39.jar \
        SE -phred33 ${raw_fastq} \
        ${trim_fastq} \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33

    fastqc -o ${fastqc_trim_dir} ${trim_fastq}
done

# run Trimmomatic and remove TruSeq adapters
adapter_input_list=("A" "J" "K" "M" "O"
        "R" "S")
for x in "${adapter_input_list[@]}"; do
    raw_fastq="${fastq_raw_dir}/input${x}.fastq.gz"
    trim_fastq="${fastq_trim_dir}/input${x}.trimmed.fastq.gz"
    java -jar ${trimmomatic_install_dir}/Trimmomatic-0.39/trimmomatic-0.39.jar \
        SE -phred33 ${raw_fastq} \
        ${trim_fastq} \
        ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33

    fastqc -o ${fastqc_trim_dir} ${trim_fastq}
done
