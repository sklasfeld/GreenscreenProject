#!/bin/bash

# 1. get down-sampled reads into bedgraph format

genome="meta/ArabidopsisGenome/TAIR10_chr_count.txt"
indir="mapped/chip/downsample"
outdir="data/bigwigs/chip/downsamp_pooled"
normalizeby=10000000 # scaling factor

mkdir -p ${outdir}

# 1A. ChIP-Seq samples

while read line; do
    replicate=$(echo $line | cut -d "," -f1)
    samp=$(echo $replicate| awk -F"_" '{rep=$NF; sub("_"rep,"",$0); print $0}')
    readSize=$(echo $line | cut -d "," -f2)
    fragSize=$(echo $line | cut -d "," -f3)
    extend=`awk -v f=${fragSize} -v r=${readSize} 'BEGIN{print f-r}'`
    orig_bam="${indir}/${replicate}.dupmark.sorted.bam"

    if [[ ! -f ${outdir}/${samp}.bigwig ]]; then
        if [[ ! -f ${outdir}/${replicate}.bg ]]; then
            # convert BAM to BED format
            bamToBed \
                -i ${orig_bam} \
                > ${outdir}/${replicate}.bed

            # sort the BED file
            sort -k 1,1 \
                ${outdir}/${replicate}.bed > \
                ${outdir}/${replicate}.sorted.bed \
                && rm ${outdir}/${replicate}.bed

            # extend the reads
            slopBed -i ${outdir}/${replicate}.sorted.bed \
                -l 0 -r ${extend} -s -g ${genome} \
                > ${outdir}/${replicate}.extend.bed \
                && rm ${outdir}/${replicate}.sorted.bed

            # normalize signal and output BEDGRAPH
            totreads=`samtools view -c ${orig_bam}`
            scaling=`awk 'BEGIN{ print '"$normalizeby"' / '"$totreads"'}'`

            genomeCoverageBed -i \
                ${outdir}/${replicate}.extend.bed \
                -g ${genome} -bg -scale $scaling | \
                awk 'BEGIN{OFS=FS="\t"} \
                {$4=sprintf("%.2f",$4)}{print}' \
                > ${outdir}/${replicate}.bg \
                && rm ${outdir}/${replicate}.extend.bed
        fi
    fi
done < meta/chip_readsize_fragsize.csv

# 1B. ChIP-Seq controls

pubControl_list=("TFL1_FD_fd_W_2020_Input_R1" "TFL1_FD_fd_W_2020_Input_R2"
        "TFL1_FD_fd_W_2020_Input_R3" "TFL1_FD_fd_W_2020_Input_R4"
        "TFL1_FD_WT_W_2020_Mock_R1" "TFL1_FD_WT_W_2020_Mock_R2"
        "TFL1_fd_W_2020_Mock_R1" "TFL1_fd_W_2020_Mock_R2"
        "LFY_W_2021_Input_R1" "LFY_W_2021_Input_R2"
        "LFY_W_2021_Input_R3" "LFY_W_2021_Mock_R1"
        "LFY_W_2021_Mock_R2" "LFY_W_2021_Mock_R3"
        "LFY_P_2016_Mock_R1" "LFY_P_2016_Mock_R2"
        "FD_S_2019_Mock_R1" "FD_ft10_tsf1_S_2019_Mock_R1"
        "TFL1_S_2020_Mock_R1" "TFL1_S_2020_Mock_R2"
        "TFL1_S_2020_Mock_R3" "LFY_P_2011_Input_R1"
        "LFY_P_2011_Input_R2" "FD_C_2020_Input_R1"
        "FD_C_2020_Input_R2")

for replicate in "${pubControl_list[@]}"; do
    samp=$(echo $replicate| awk -F"_" '{rep=$NF; sub("_"rep,"",$0); print $0}')
    if [[ ! -f ${outdir}/${samp}.bigwig ]]; then
        if [[ ! -f ${outdir}/${replicate}.bg ]]; then
            orig_bam="${indir}/${replicate}.dupmark.sorted.bam"

            # normalize signal and output BEDGRAPH
            totreads=`samtools view -c ${orig_bam}`
            scaling=`awk 'BEGIN{ print '"$normalizeby"' / '"$totreads"'}'`

            genomeCoverageBed -ibam ${orig_bam} \
                -bg -scale $scaling | \
                awk 'BEGIN{OFS=FS="\t"} \
                {$4=sprintf("%.2f",$4)}{print}' \
                > ${outdir}/${replicate}.bg
        fi
    fi
done

# 2. Get the average signal of the replicates

# 2A. return average signal of one rep
pool_one=("FD_ft10_tsf1_S_2019_Mock" "FD_S_2019_Mock")
for samp in "${pool_one[@]}"; do
    if [[ ! -f ${outdir}/${samp}.bigwig ]]; then
        cp ${outdir}/${samp}_R1.bg \
        ${outdir}/${samp}.bigwig
    fi
done

# 2B. return average signal of two reps
pool_two=("FD_W_2020" "TFL1_fd_W_2020"
  "LFY_P_2016" "FD_S_2019"
  "FD_ft10_tsf1_S_2019" "LFY_P_2011"
  "FD_C_2020" "FD_C_2020_Input"
  "LFY_P_2011_Input" "LFY_P_2016_Mock"
  "TFL1_fd_W_2020_Mock" "TFL1_FD_WT_W_2020_Mock")
for samp in "${pool_two[@]}"; do
    if [[ ! -f ${outdir}/${samp}.bigwig ]]; then
        bedtools unionbedg -i \
        ${outdir}/${samp}_R1.bg \
        ${outdir}/${samp}_R2.bg \
        > ${outdir}/${samp}.unionbg \
        && rm ${outdir}/${samp}_R1.bg \
        ${outdir}/${samp}_R2.bg

        awk 'BEGIN{OFS="\t"} \
        $1!="ChrC" && $1!="ChrM"{ \
        avg=($4+$5)/2;
        print $1,$2,$3,avg}' \
        ${outdir}/${samp}.unionbg > \
        ${outdir}/${samp}.bedgraph \
        && rm ${outdir}/${samp}.unionbg

        bedGraphToBigWig ${outdir}/${samp}.bedgraph \
        ${genome} ${outdir}/${samp}.bigwig \
        && rm ${outdir}/${samp}.bedgraph
    fi
done

# 2C. return average signal of three reps
pool_three=("TFL1_S_2020" "LFY_W_2021"
    "LFY_W_2021_Input" "LFY_W_2021_Mock"
    "TFL1_S_2020_Mock")
for samp in "${pool_three[@]}"; do
    if [[ ! -f ${outdir}/${samp}.bigwig ]]; then
        bedtools unionbedg -i \
            ${outdir}/${samp}_R1.bg \
            ${outdir}/${samp}_R2.bg \
            ${outdir}/${samp}_R3.bg \
            > ${outdir}/${samp}.unionbg \
            && rm ${outdir}/${samp}_R1.bg \
            ${outdir}/${samp}_R2.bg \
            ${outdir}/${samp}_R3.bg

        awk 'BEGIN{OFS="\t"} \
            $1!="ChrC" && $1!="ChrM"{ \
            avg=($4+$5+$6)/3;
            print $1,$2,$3,avg}' \
            ${outdir}/${samp}.unionbg > \
            ${outdir}/${samp}.bedgraph \
            && rm ${outdir}/${samp}.unionbg

        bedGraphToBigWig ${outdir}/${samp}.bedgraph \
            ${genome} ${outdir}/${samp}.bigwig \
            && rm ${outdir}/${samp}.bedgraph
    fi
done

# 2D. return average signal of four reps
pool_four=("TFL1_FD_fd_W_2020_Input" "TFL1_W_2020")
for samp in "${pool_four[@]}"; do
    if [[ ! -f ${outdir}/${samp}.bigwig ]]; then
        bedtools unionbedg -i \
            ${outdir}/${samp}_R1.bg \
            ${outdir}/${samp}_R2.bg \
            ${outdir}/${samp}_R3.bg \
            ${outdir}/${samp}_R4.bg \
            > ${outdir}/${samp}.unionbg \
            && rm ${outdir}/${samp}_R1.bg \
            ${outdir}/${samp}_R2.bg \
            ${outdir}/${samp}_R3.bg \
            ${outdir}/${samp}_R4.bg

        awk 'BEGIN{OFS="\t"} \
            $1!="ChrC" && $1!="ChrM"{ \
            avg=($4+$5+$6+$7)/4;
            print $1,$2,$3,avg}' \
            ${outdir}/${samp}.unionbg > \
            ${outdir}/${samp}.bedgraph \
            && rm ${outdir}/${samp}.unionbg

        bedGraphToBigWig ${outdir}/${samp}.bedgraph \
            ${genome} ${outdir}/${samp}.bigwig \
            && rm ${outdir}/${samp}.bedgraph
    fi
done