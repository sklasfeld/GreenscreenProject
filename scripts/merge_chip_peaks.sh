#!/bin/bash

chip_peaks_dir="data/macs2_out/chipPeaks/gsMask_qval10"

cat \
	${chip_peaks_dir}/FD_W_2020_R1_peaks.broadPeak \
	${chip_peaks_dir}/FD_W_2020_R2_peaks.broadPeak \
	${chip_peaks_dir}/TFL1_A_W_2020_R1_peaks.broadPeak \
	${chip_peaks_dir}/TFL1_A_W_2020_R2_peaks.broadPeak \
	${chip_peaks_dir}/TFL1_B_W_2020_R1_peaks.broadPeak \
	${chip_peaks_dir}/TFL1_B_W_2020_R2_peaks.broadPeak \
	${chip_peaks_dir}/TFL1_fd_W_2020_R1_peaks.broadPeak \
	${chip_peaks_dir}/TFL1_fd_W_2020_R2_peaks.broadPeak \
	${chip_peaks_dir}/LFY_W_2021_R1_peaks.broadPeak \
	${chip_peaks_dir}/LFY_W_2021_R2_peaks.broadPeak \
	${chip_peaks_dir}/LFY_W_2021_R3_peaks.broadPeak \
	${chip_peaks_dir}/LFY_P_2016_R1_peaks.broadPeak \
	${chip_peaks_dir}/LFY_P_2016_R2_peaks.broadPeak \
	${chip_peaks_dir}/FD_S_2019_R1_peaks.broadPeak \
	${chip_peaks_dir}/FD_S_2019_R2_peaks.broadPeak \
	${chip_peaks_dir}/FD_ft10_tsf1_S_2019_R1_peaks.broadPeak \
	${chip_peaks_dir}/FD_ft10_tsf1_S_2019_R2_peaks.broadPeak \
	${chip_peaks_dir}/TFL1_S_2020_R1_peaks.broadPeak \
	${chip_peaks_dir}/TFL1_S_2020_R2_peaks.broadPeak \
	${chip_peaks_dir}/TFL1_S_2020_R3_peaks.broadPeak \
	${chip_peaks_dir}/LFY_P_2011_R1_peaks.broadPeak \
	${chip_peaks_dir}/LFY_P_2011_R2_peaks.broadPeak \
	${chip_peaks_dir}/FD_C_2020_R1_peaks.broadPeak \
	${chip_peaks_dir}/FD_C_2020_R2_peaks.broadPeak | \
	sort -k1,1 -k2,2n | \
	bedtools merge -i - > \
	${chip_peaks_dir}/ChIPseq_Peaks.merged.bed
