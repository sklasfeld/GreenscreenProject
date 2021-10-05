#!/bin/bash

in_dir="mapped/chip"
out_dir="mapped/chip/downsample"
JVARKIT_PATH="/usr/src/jvarkit"
mkdir -p ${out_dir}

# There is no need to down-sample anything
# with a single replicate
single_rep=("FD_S_2019_Mock" 
"FD_ft10_tsf1_S_2019_Mock")
for samp in "${single_rep[@]}"; do
  cp ${in_dir}/${samp}_R1.dupmark.sorted.bam \
    ${out_dir}/${samp}_R1.dupmark.sorted.bam
  cp ${in_dir}/${samp}_R1.dupmark.sorted.bam.bai \
    ${out_dir}/${samp}_R1.dupmark.sorted.bam.bai
done

# down-sample given two reps
pool_two=("FD_W_2020" "TFL1_fd_W_2020"
  "LFY_P_2016" "FD_S_2019" 
  "FD_ft10_tsf1_S_2019" "LFY_P_2011"
  "FD_C_2020" "TFL1_FD_WT_W_2020_Mock" 
  "TFL1_fd_W_2020_Mock"
  "LFY_P_2016_Mock" "LFY_P_2011_Input"
  "FD_C_2020_Input")
for samp in "${pool_two[@]}"; do
  depth1=`samtools view -c \
    ${in_dir}/${samp}_R1.dupmark.sorted.bam`
  depth2=`samtools view -c \
    ${in_dir}/${samp}_R2.dupmark.sorted.bam`
  arrName=(${depth1} ${depth2})
  minIndex "${arrName[@]}"
  let "min_idx = $min_idx + 1"
  for (( rep=1; rep<=2; rep++ )); do
    if [[ ${min_idx} -eq ${rep} ]]; then
        # copy this replicate which
        # has the smallest read depth
        cp ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam \
          ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam
        cp ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam.bai \
          ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam.bai
    else
        # downsample these replicates
        java -jar ${JVARKIT_PATH}/dist/biostar145820.jar \
            -n ${min_val} \
            -o ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam \
            ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam
    fi
  done
done

# down-sample given three reps
pool_three=("TFL1_S_2020" "LFY_W_2021"
"LFY_W_2021_Input" "LFY_W_2021_Mock"
"TFL1_S_2020_Mock")
for samp in "${pool_three[@]}"; do
  depth1=`samtools view -c \
    ${in_dir}/${samp}_R1.dupmark.sorted.bam`
  depth2=`samtools view -c \
    ${in_dir}/${samp}_R2.dupmark.sorted.bam`
  depth3=`samtools view -c \
    ${in_dir}/${samp}_R3.dupmark.sorted.bam`
  arrName=(${depth1} ${depth2} ${depth3})
  minIndex "${arrName[@]}"
  let "min_idx = $min_idx + 1"
  for (( rep=1; rep<=3; rep++ )); do
    if [[ ${min_idx} -eq ${rep} ]]; then
        # copy this replicate which
        # has the smallest read depth
        cp ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam \
          ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam
        cp ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam.bai \
          ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam.bai
    else
        # downsample these replicates
        java -jar ${JVARKIT_PATH}/dist/biostar145820.jar \
            -n ${min_val} \
            -o ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam \
            ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam
    fi
  done
done

# down-sample given four reps
pool_four=("TFL1_W_2020" "TFL1_FD_fd_W_2020_Input")
for samp in "${pool_four[@]}"; do
  depth1=`samtools view -c \
    ${in_dir}/${samp}_R1.dupmark.sorted.bam`
  depth2=`samtools view -c \
    ${in_dir}/${samp}_R2.dupmark.sorted.bam`
  depth3=`samtools view -c \
    ${in_dir}/${samp}_R3.dupmark.sorted.bam`
  depth4=`samtools view -c \
    ${in_dir}/${samp}_R4.dupmark.sorted.bam`
  arrName=(${depth1} ${depth2} ${depth3} ${depth4})
  minIndex "${arrName[@]}"
  let "min_idx = $min_idx + 1"
  for (( rep=1; rep<=4; rep++ )); do
    if [[ ${min_idx} -eq ${rep} ]]; then
        # copy this replicate which
        # has the smallest read depth
        cp ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam \
          ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam
    else
        # downsample these replicates
        java -jar ${JVARKIT_PATH}/dist/biostar145820.jar \
            -n ${min_val} \
            -o ${out_dir}/${samp}_R${rep}.dupmark.sorted.bam \
            ${in_dir}/${samp}_R${rep}.dupmark.sorted.bam
    fi
  done
done
