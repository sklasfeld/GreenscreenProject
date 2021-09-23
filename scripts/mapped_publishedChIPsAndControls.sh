#!/bin/bash

nthreads=2
PICARD_PATH="/usr/src/picard/build/libs"

# create output directory
mkdir -p mapped/chip

while read samp; do
    # run bowtie2
    bowtie2  --phred33 -q \
	-x meta/ArabidopsisGenome/bowtie2_genome_dir/TAIR10 \
        -S mapped/chip/${samp}.sam \
        fastq/trimmed/${samp}.trimmed.fastq.gz
    # sort the reads
    samtools sort -o mapped/chip/${samp}.bam \
        mapped/chip/${samp}.sam
    # index the bam file
    samtools index mapped/chip/${samp}.bam \
    # remove reads without MAPQ>=30
    samtools view -@ ${nthreads} -F 772 -q 30 \
        -b mapped/chip/${samp}.bam \
	Chr1 Chr2 Chr3 Chr4 Chr5 | \
        samtools sort - -o mapped/chip/${samp}.filter.bam
    # index filtered reads
    samtools index mapped/chip/${samp}.filter.bam
    # mark duplicates with picard
    java -jar ${PICARD_PATH}/picard.jar MarkDuplicates \
        I=mapped/chip/${samp}.filter.bam \
        O=mapped/chip/${samp}.dupmark.bam \
        M=mapped/chip/${samp}.dup.qc \ 
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=false ASSUME_SORTED=true
    # sort reads after marking the duplicates
    samtools sort -o mapped/chip/${samp}.dupmark.sorted.bam \
        mapped/chip/${samp}.dupmark.bam
    # index the sorted reads
    samtools index mapped/chip/${samp}.dupmark.sorted.bam
    # remove duplicates
    samtools view -@ ${nthreads} -F 1796 -q 30 \
        -b mapped/chip/${samp}.dupmark.sorted.bam | \
        samtools sort - -o mapped/chip/${samp}.noDups.bam
    # index unique reads
    samtools index mapped/chip/${samp}.noDups.bam
done < meta/chip_controls_sampleIDs.list
