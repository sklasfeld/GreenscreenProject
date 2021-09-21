#!/bin/bash

genome="meta/ArabidopsisGenome/TAIR10_chr_count.txt"
working_dir="mapped/input"
normalizeby=10000000 # scaling factor

input_list=("A" "B" "C" "D" "E"
        "F" "G" "H" "I" "J" "K"
        "L" "M" "N" "O" "P" "Q"
        "R" "S" "T")

for x in "${input_list[@]}"; do
    orig_bam="${working_dir}/input${x}.dupmark.sorted.bam"

    # normalize signal and output BEDGRAPH
    totreads=`samtools view -c ${orig_bam}`
    scaling=`awk 'BEGIN{ print '"$normalizeby"' / '"$totreads"'}'`

    genomeCoverageBed -ibam ${orig_bam} \
      -bg -scale $scaling | \
      awk 'BEGIN{OFS=FS="\t"} \
      {$4=sprintf("%.2f",$4)}{print}' \
      > ${working_dir}/${samp}.bg

    # compress BEDGRAPH to BIGWIG FORMAT
    bedGraphToBigWig ${working_dir}/${samp}.bg \
      ${genome} ${working_dir}/${samp}.bw \
    && rm ${working_dir}/${samp}.bg
done
