#!/bin/bash

if [ "$#" -ne 3 ]; then
        echo -e "\n> Goal: Calculate genome contig sizes"
        echo -e "\n> Command Format:"
        echo -e "\nbash /path/to/generate_20input_greenscreenBed.sh [qval] [merge_distance] [distinct_ninputs]"
        echo -e "\n> Parameters:"
        echo -e "\n\t* qval: average basepair q-value threshold (log10) (eg. set to 10 for q<=10^-10)\n"
        echo -e "\n\t* merge_distance: neighboring regions within this distance are merged into one\n"
        echo -e "\n\t* distinct_ninputs: final regions must also show signal in at least this number of inputs\n"
else

	qval=$1
	merge_distance=$2
	distinct_ninputs=$3

	# directory to output MACS2 results
	macs_out="data/macs2_out/inputControls"

	# make output directory
	out_dir="${macs_out}/qval${qval}"
	mkdir -p ${out_dir}

	input_list=("A" "B" "C" "D" "E"
        	"F" "G" "H" "I" "J" "K"
        	"L" "M" "N" "O" "P" "Q"
	        "R" "S" "T")

	for x in "${input_list[@]}"; do
		# remove all peaks that do not have an
		# average base pair q-value <=10^(-${q})
    		awk -F"\t" -v q=${qval} 'BEGIN{OFS="\t"} \
        		$9>=q && $1!="ChrC" && $1!="ChrM"{print}' \
        		${macs_out}/input${x}_peaks.broadPeak > \
        		${out_dir}/input${x}_peaks.broadPeak
	done

	concat_file="${out_dir}/concat_input_peaks.broadPeak"
	merge_file="${out_dir}/merge${merge_distance}bp_20inputs.txt"
	final_greenscreen="${out_dir}/gs_merge${merge_distance}bp_call${distinct_ninputs}_20inputs.txt"

	# concatenate peaks and set column 4 (peak name) to the sample ID
	cat ${out_dir}/inputA_peaks.broadPeak \
		${out_dir}/inputB_peaks.broadPeak \
		${out_dir}/inputC_peaks.broadPeak \
		${out_dir}/inputD_peaks.broadPeak \
		${out_dir}/inputE_peaks.broadPeak \
		${out_dir}/inputF_peaks.broadPeak \
		${out_dir}/inputG_peaks.broadPeak \
		${out_dir}/inputH_peaks.broadPeak \
		${out_dir}/inputI_peaks.broadPeak \
		${out_dir}/inputJ_peaks.broadPeak \
		${out_dir}/inputK_peaks.broadPeak \
		${out_dir}/inputL_peaks.broadPeak \
		${out_dir}/inputM_peaks.broadPeak \
		${out_dir}/inputN_peaks.broadPeak \
		${out_dir}/inputO_peaks.broadPeak \
		${out_dir}/inputP_peaks.broadPeak \
		${out_dir}/inputQ_peaks.broadPeak \
		${out_dir}/inputR_peaks.broadPeak \
		${out_dir}/inputS_peaks.broadPeak \
		${out_dir}/inputT_peaks.broadPeak | \
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
	    ${merge_file} \
	    > ${final_greenscreen}
fi
