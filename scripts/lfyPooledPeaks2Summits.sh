#!/bin/bash

in_dir="data/macs2_out/chipPeaks/gsMask_qval10"

# generate summit BED file
$ awk -F"\t" 'BEGIN{OFS="\t"} \
    {print $1,$2+$10,$2+$10+1,$4,$5,$6}' \
    ${in_dir}/LFY_W_2021_peaks.narrowPeak > \
    ${in_dir}/LFY_W_2021_summits.bed


# generate summit NARROWPEAK file
$ awk -F"\t" 'BEGIN{OFS="\t"} \
    {print $1,$2+$10,$2+$10+1,$4,$5,$6,$7,$8,$9,"0"}' \
    ${in_dir}/LFY_W_2021_peaks.narrowPeak > \
    ${in_dir}/LFY_W_2021_summits.narrowPeak
