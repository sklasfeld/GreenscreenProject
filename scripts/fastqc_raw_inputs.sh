#!/bin/bash

input_list=("A" "B" "C" "D" "E"
        "F" "G" "H" "I" "J" "K"
        "L" "M" "N" "O" "P" "Q"
        "R" "S" "T")

for x in "${input_list[@]}"; do
    fastqc -o data/FASTQC/raw fastq/raw/input${x}.fastq
done
