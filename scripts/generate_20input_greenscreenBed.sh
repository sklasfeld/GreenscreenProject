#!/bin/bash

merge_distance=1000
distinct_ninputs=10

concat_file="data/macs2_out/inputControls/qval10/concat_input_peaks.broadPeak"
merge_file="data/macs2_out/inputControls/qval10/merged_input_peaks.txt"
final_greenscreen="data/arabidopsis_greenscreen_20inputs.bed"

# concatenate peaks and set column 4 (peak name) to the sample ID
cat data/macs2_out/inputControls/qval10/inputA_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputB_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputC_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputD_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputE_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputF_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputG_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputH_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputI_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputJ_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputK_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputL_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputM_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputN_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputO_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputP_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputQ_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputR_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputS_peaks.broadPeak \
	data/macs2_out/inputControls/qval10/inputT_peaks.broadPeak | \
	awk 'BEGIN{OFS="\t"} \
	{sub(/(_peak_)[0-9]+/,"",$4); print $0}' | \
	sort -k1,1 -k2,2n \
	> ${concat_file}

# merge all overlapping regions and
# regions within ${merge_distance} bp apart
bedtools merge -c 4,5,6,7,8,9,4 \
    -o distinct,max,distinct,max,max,max,count_distinct \
	-i ${concat_file} -d ${merge_distance} > ${merge_file}

# filter out regions that are called in less
# than ${distinct_ninputs} distinct samples
awk -v thresh=${distinct_ninputs} -F "\t" \
    'BEGIN{OFS="\t"} $10>=thresh{print}' \
    > ${final_greenscreen}
