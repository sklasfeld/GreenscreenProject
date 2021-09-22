#!/bin/bash

in_dir="mapped/chip"
out_dir="mapped/chip/downsample"

mkdir -p ${out_dir}

# seed (this is set for reproducibility purposes)
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

# function to down-sample reads
downsamp (){
  idir="$1"
  odir="$2"
  sampl="$3"
  n="$4"

  # get header for output sam file
  samtools view -H \
    ${idir}/${sampl}.dupmark.sorted.bam \
    > ${odir}/${sampl}.dupmark.sorted.sam
  # 1. convert bam file to human readable sam format
  # 2. shuffle the reads
  # 3. return the first  $n reads
  samtools view \
    ${idir}/${sampl}.dupmark.sorted.bam | \
    shuf \
    --random-source=<(get_seeded_random 42) - | \
    head -n ${n} \
    >> ${odir}/${sampl}.dupmark.sorted.sam

  # convert the sam file to bam format
  samtools sort -O BAM \
    ${odir}/${sampl}.dupmark.sorted.sam \
    > ${odir}/${sampl}.dupmark.sorted.bam
  samtools index \
    ${odir}/${sampl}.dupmark.sorted.bam
}

# function to find the minimum
# value in an array
minIndex(){
   arr=("$@")
   min_val=${arr[0]}
   min_idx=0
   for i in ${!arr[@]}; do
        cur_val=${arr[${i}]}
        if [[ ${cur_val} -lt ${min_val} ]]; then
                min_val=${arr[$i]}
                min_idx=${i}
        fi
   done

}

# There is no need to down-sample anything
# with a single replicate
single_rep=("FD_S_2019_Mock" 
"FD_ft10_tsf1_S_2019_Mock")
for samp in "${single_rep[@]}"; do
  cp ${in_dir}/${samp}_R1.dupmark.sorted.bam \
    ${out_dir}/${samp}_R1.dupmark.sorted.bam
done

# down-sample given two reps
pool_two=("FD_W_2020" "TFL1_A_W_2020"
  "TFL1_B_W_2020" "TFL1_fd_W_2020"
  "LFY_P_2016" "FD_S_2019" 
  "FD_ft10_tsf1_S_2019" "LFY_P_2011"
  "FD_C_2020" "TFL1_FD_W_2020_Mock" 
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
    else
      # downsample these replicates
      downsamp ${in_dir} ${out_dir} ${samp}_R${rep} ${min_val}
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
    else
      # downsample these replicates
      downsamp ${in_dir} ${out_dir} ${samp}_R${rep} ${min_val}
    fi
  done
done

# down-sample given four reps
pool_four=("TFL1_FD_fd_W_2020_Input")
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
      downsamp ${in_dir} ${out_dir} ${samp}_R${rep} ${min_val}
    fi
  done
done
