#!/bin/bash

genome="meta/ArabidopsisGenome/TAIR10_chr_count.txt"
working_dir="mapped/chip"
normalizeby=10000000 # scaling factor

while read line; do
    samp=$(echo $line | cut -d "," -f1)
    extend=$(echo $line | cut -d "," -f2)
    orig_bam="${working_dir}/${samp}.dupmark.sorted.bam"
    # convert BAM to BED format
    bamToBed \
      -i ${orig_bam} \
      > ${working_dir}/${samp}.bed

    # sort the BED file
    sort -k 1,1 \
      ${working_dir}/${samp}.bed > \
      ${working_dir}/${samp}.sorted.bed

    # extend the reads
    slopBed -i ${working_dir}/${samp}.sorted.bed \
      -l 0 -r ${extend} -s -g ${genome} \
      > ${working_dir}/${samp}.extend.bed

    # normalize signal and output BEDGRAPH
    totreads=`samtools view -c ${orig_bam}`
    scaling=`awk 'BEGIN{ print '"$normalizeby"' / '"$totreads"'}'`
    
    genomeCoverageBed -i \
      ${working_dir}/${samp}.extend.bed \
      -g ${genome} -bg -scale $scaling | \
      awk 'BEGIN{OFS=FS="\t"} \
      {$4=sprintf("%.2f",$4)}{print}' \
      > ${working_dir}/${samp}.bg

    # compress BEDGRAPH to BIGWIG FORMAT
    bedGraphToBigWig ${working_dir}/${samp}.bg \
      ${genome} ${working_dir}/${samp}.bw
done < meta/chipFragSize.csv
