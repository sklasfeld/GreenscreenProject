#!/bin/bash

working_dir="mapped/chip"
genome="meta/ArabidopsisGenome/TAIR10_chr_count.txt"

# return average signal of two reps
pool_two=("FD_W_2020" "TFL1_A_W_2020"
  "TFL1_B_W_2020" "TFL1_fd_W_2020"
  "LFY_P_2016" "FD_S_2019" 
  "FD_ft10_tsf1_S_2019" "LFY_P_2011"
  "FD_C_2020")
for samp in "${pool_two[@]}"; do
  bedtools unionbedg -i \
    ${working_dir}/${samp}_R1.bg \
    ${working_dir}/${samp}_R2.bg \
    > ${working_dir}/${samp}.unionbg
  awk 'BEGIN{OFS="\t"} \
    $1!="ChrC" && $1!="ChrM"{ \
    avg=($4+$5)/2;
    print $1,$2,$3,avg}' \
    ${working_dir}/${samp}.unionbg > \
    ${working_dir}/${samp}.bedgraph
  bedGraphToBigWig ${working_dir}/${samp}.bedgraph \
    ${genome} ${working_dir}/${samp}.bigwig
done

# return average signal of three reps
pool_three=("TFL1_S_2020" "LFY_W_2021")
for samp in "${pool_three[@]}"; do
  bedtools unionbedg -i \
    ${working_dir}/${samp}_R1.bg \
    ${working_dir}/${samp}_R2.bg \
    ${working_dir}/${samp}_R3.bg \
    > ${working_dir}/${samp}.unionbg
  awk 'BEGIN{OFS="\t"} \
    $1!="ChrC" && $1!="ChrM"{ \
    avg=($4+$5+$6)/3;
    print $1,$2,$3,avg}' \
    ${working_dir}/${samp}.unionbg > \
    ${working_dir}/${samp}.bedgraph
  bedGraphToBigWig ${working_dir}/${samp}.bedgraph \
    ${genome} ${working_dir}/${samp}.bigwig
done
