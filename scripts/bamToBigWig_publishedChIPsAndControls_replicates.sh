#!/bin/bash

genome="meta/ArabidopsisGenome/TAIR10_chr_count.txt"
indir="mapped/chip"
outdir="data/bigwigs/chip/individual_replicates"
normalizeby=10000000 # scaling factor

mkdir -p ${outdir}

# ChIP-Seq samples

while read line; do
    samp=$(echo $line | cut -d "," -f1)
    readSize=$(echo $line | cut -d "," -f2)
    fragSize=$(echo $line | cut -d "," -f3)
    extend=`awk -v f=${fragSize} -v r=${readSize} 'BEGIN{print f-r}'`
    orig_bam="${indir}/${samp}.dupmark.sorted.bam"
    if [[ ! -f "${outdir}/${samp}.bw" ]]; then
        # convert BAM to BED format
        bamToBed \
          -i ${orig_bam} \
          > ${outdir}/${samp}.bed

        # sort the BED file
        sort -k 1,1 \
          ${outdir}/${samp}.bed > \
          ${outdir}/${samp}.sorted.bed

        # extend the reads
        slopBed -i ${outdir}/${samp}.sorted.bed \
            -l 0 -r ${extend} -s -g ${genome} \
            > ${outdir}/${samp}.extend.bed

        # normalize signal and output BEDGRAPH
        totreads=`samtools view -c ${orig_bam}`
        scaling=`awk 'BEGIN{ print '"$normalizeby"' / '"$totreads"'}'`

        genomeCoverageBed -i \
           ${outdir}/${samp}.extend.bed \
           -g ${genome} -bg -scale $scaling | \
           awk 'BEGIN{OFS=FS="\t"} \
           {$4=sprintf("%.2f",$4)}{print}' \
           > ${outdir}/${samp}.bg

        # compress BEDGRAPH to BIGWIG FORMAT
        bedGraphToBigWig ${outdir}/${samp}.bg \
           ${genome} ${outdir}/${samp}.bw
    fi
done < meta/chip_readsize_fragsize.csv

# published ChIP-Seq controls

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

for samp in "${pubControl_list[@]}"; do
    if [[ ! -f "${outdir}/${samp}.bw" ]]; then
        orig_bam="${indir}/${samp}.dupmark.sorted.bam"

        # normalize signal and output BEDGRAPH
        totreads=`samtools view -c ${orig_bam}`
        scaling=`awk 'BEGIN{ print '"$normalizeby"' / '"$totreads"'}'`

        genomeCoverageBed -ibam ${orig_bam} \
          -bg -scale $scaling | \
          awk 'BEGIN{OFS=FS="\t"} \
          {$4=sprintf("%.2f",$4)}{print}' \
          > ${outdir}/${samp}.bg

        # compress BEDGRAPH to BIGWIG FORMAT
        bedGraphToBigWig ${outdir}/${samp}.bg \
          ${genome} ${outdir}/${samp}.bw \
          && rm ${outdir}/${samp}.bg
    fi
done
