#!/bin/bash

# 1. get down-sampled reads into bedgraph format

genome="meta/ArabidopsisGenome/TAIR10_chr_count.txt"
working_dir="mapped/chip/downsample"
normalizeby=10000000 # scaling factor

while read line; do
    replicate=$(echo $line | cut -d "," -f1)
    samp=$(echo $replicate| awk -F"_" '{rep=$NF; sub("_"rep,"",$0); print $0}')
    readSize=$(echo $line | cut -d "," -f2)
    fragSize=$(echo $line | cut -d "," -f3)
    extend=`awk -v f=${fragSize} -v r=${readSize} 'BEGIN{print f-r}'`
    orig_bam="${working_dir}/${replicate}.dupmark.sorted.bam"
    
    if [[ ! -f ${working_dir}/${samp}.bigwig ]]; then
        # convert BAM to BED format
        bamToBed \
            -i ${orig_bam} \
            > ${working_dir}/${replicate}.bed

        # sort the BED file
        sort -k 1,1 \
            ${working_dir}/${replicate}.bed > \
            ${working_dir}/${replicate}.sorted.bed \
            && rm ${working_dir}/${replicate}.bed

        # extend the reads
        slopBed -i ${working_dir}/${replicate}.sorted.bed \
            -l 0 -r ${extend} -s -g ${genome} \
            > ${working_dir}/${replicate}.extend.bed \
            && rm ${working_dir}/${replicate}.sorted.bed

        # normalize signal and output BEDGRAPH
        totreads=`samtools view -c ${orig_bam}`
        scaling=`awk 'BEGIN{ print '"$normalizeby"' / '"$totreads"'}'`

        genomeCoverageBed -i \
            ${working_dir}/${replicate}.extend.bed \
            -g ${genome} -bg -scale $scaling | \
            awk 'BEGIN{OFS=FS="\t"} \
            {$4=sprintf("%.2f",$4)}{print}' \
            > ${working_dir}/${replicate}.bg \
            && rm ${working_dir}/${replicate}.extend.bed
    fi
done < meta/chip_readsize_fragsize.csv


# 2. return average signal of two reps
pool_two=("FD_W_2020" "TFL1_A_W_2020"
  "TFL1_B_W_2020" "TFL1_fd_W_2020"
  "LFY_P_2016" "FD_S_2019"
  "FD_ft10_tsf1_S_2019" "LFY_P_2011"
  "FD_C_2020")
for samp in "${pool_two[@]}"; do
    if [[ ! -f ${working_dir}/${samp}.bigwig ]]; then
        bedtools unionbedg -i \
        ${working_dir}/${samp}_R1.bg \
        ${working_dir}/${samp}_R2.bg \
        > ${working_dir}/${samp}.unionbg \
        && rm ${working_dir}/${samp}_R1.bg \
        ${working_dir}/${samp}_R2.bg

        awk 'BEGIN{OFS="\t"} \
        $1!="ChrC" && $1!="ChrM"{ \
        avg=($4+$5)/2;
        print $1,$2,$3,avg}' \
        ${working_dir}/${samp}.unionbg > \
        ${working_dir}/${samp}.bedgraph \
        && rm ${working_dir}/${samp}.unionbg

        bedGraphToBigWig ${working_dir}/${samp}.bedgraph \
        ${genome} ${working_dir}/${samp}.bigwig \
        && rm ${working_dir}/${samp}.bedgraph
    fi
done

# 3. return average signal of three reps
pool_three=("TFL1_S_2020" "LFY_W_2021")
for samp in "${pool_three[@]}"; do
    if [[ ! -f ${working_dir}/${samp}.bigwig ]]; then
        bedtools unionbedg -i \
            ${working_dir}/${samp}_R1.bg \
            ${working_dir}/${samp}_R2.bg \
            ${working_dir}/${samp}_R3.bg \
            > ${working_dir}/${samp}.unionbg \
            && rm ${working_dir}/${samp}_R1.bg \
            ${working_dir}/${samp}_R2.bg \
            ${working_dir}/${samp}_R3.bg

        awk 'BEGIN{OFS="\t"} \
            $1!="ChrC" && $1!="ChrM"{ \
            avg=($4+$5+$6)/3;
            print $1,$2,$3,avg}' \
            ${working_dir}/${samp}.unionbg > \
            ${working_dir}/${samp}.bedgraph \
            && rm ${working_dir}/${samp}.unionbg

        bedGraphToBigWig ${working_dir}/${samp}.bedgraph \
            ${genome} ${working_dir}/${samp}.bigwig \
            && rm ${working_dir}/${samp}.bedgraph
    fi
done