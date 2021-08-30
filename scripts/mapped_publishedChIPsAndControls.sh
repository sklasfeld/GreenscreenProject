#!/bin/bash

nthreads=2
while read samp; do
    # run bowtie2
    bowtie2 -x meta/Araport11/bowtie2_genome_dir/TAIR10 \
        -S mapped/chip/${samp}.sam \
        fastq/trimmed/${samp}.trimmed.fastq
    # sort the reads
    samtools sort -o mapped/chip/${samp}.bam \
        mapped/chip/${samp}.sam
    # index the bam file
    samtools index mapped/chip/${samp}.bam \
    # remove reads without MAPQ>=30
    samtools view -@ ${nthreads} -F 1796 -q 30 \
        -b mapped/chip/${samp}.bam | \
        samtools sort - -o mapped/chip/${samp}.filter.bam
    # index filtered reads
    samtools index mapped/chip/${samp}.filter.bam
    # mark duplicates with picard
    java -jar picard.jar MarkDuplicates \
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
