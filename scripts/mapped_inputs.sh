#!/bin/bash

nthreads=2
PICARD_PATH="/usr/src/picard/build/libs"
input_list=("A" "B" "C" "D" "E"
        "F" "G" "H" "I" "J" "K"
        "L" "M" "N" "O" "P" "Q"
        "R" "S" "T")

for x in "${input_list[@]}"; do
    # run bowtie2
    bowtie2 -x meta/ArabidopsisGenome/bowtie2_genome_dir/TAIR10 \
        -S mapped/input/input${x}.sam \
        fastq/trimmed/input${x}.trimmed.fastq
    # sort the reads
    samtools sort -o mapped/input/input${x}.bam \
        mapped/input/input${x}.sam
    # index the bam file
    samtools index mapped/input/input${x}.bam
    # remove reads without MAPQ>=30
    samtools view -@ ${nthreads} -F 1796 -q 30 \
        -b mapped/input/input${x}.bam | \
        samtools sort - -o mapped/input/input${x}.filter.bam
    # index filtered reads
    samtools index mapped/input/input${x}.filter.bam
    # mark duplicates with picard
    java -jar ${PICARD_PATH}/picard.jar MarkDuplicates \
        I=mapped/input/input${x}.filter.bam \
        O=mapped/input/input${x}.dupmark.bam \
        M=mapped/input/input${x}.dup.qc \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=false ASSUME_SORTED=true
    # sort reads after marking the duplicates
    samtools sort -o mapped/input/input${x}.dupmark.sorted.bam \
        mapped/input/input${x}.dupmark.bam
    # index the sorted reads
    samtools index mapped/input/input${x}.dupmark.sorted.bam
done
