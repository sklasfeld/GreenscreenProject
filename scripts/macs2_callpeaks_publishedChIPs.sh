#!/bin/bash

raw_bam_dir="mapped/chip"
downsamp_bam_dir="${raw_bam_dir}/downsample"
bam_suffix="dupmark.sorted.bam"
macs2_out="data/macs2_out/chipPeaks"
gs_regions="data/macs2_out/inputControls/qval10/gs_merge5000bp_call10_20inputs.txt"
# average basepair q-value threshold (log10)
q=10

# make macs2 output directory
mkdir -p ${macs2_out}/noMask_qval${q}
mkdir -p ${macs2_out}/gsMask_qval${q}

while read line; do
    chip_name=$(echo $line | cut -d "," -f1)
    chip_nreps=$(echo $line | cut -d "," -f2)
    chip_fsize=$(echo $line | cut -d "," -f3)
    cntl_name=$(echo $line | cut -d "," -f4)
    cntl_nreps=$(echo $line | cut -d "," -f5)

    c_param="${downsamp_bam_dir}/${cntl_name}_R1.${bam_suffix}"
    for ((c=2; c<=${cntl_nreps}; c++ )); do
        c_param="${c_param} ${downsamp_bam_dir}/${cntl_name}_R${c}.${bam_suffix}"
    done

    # call peaks on individual replicates

    # run MACS2
    for ((t=1; t<=${chip_nreps}; t++ )); do
        t_raw_rep="${raw_bam_dir}/${chip_name}_R${t}.${bam_suffix}"

        final_peaks="${macs2_out}/gsMask_qval${q}/${chip_name}_R${t}_peaks.narrowPeak"
        if [[ ! -f $final_peaks ]]; then
            macs2 callpeak -t ${t_raw_rep} \
                -c ${c_param} \
                -f BAM --keep-dup auto \
                --nomodel --extsize ${chip_fsize} -g 101274395 \
                --outdir ${macs2_out} -n ${chip_name}_R${t}


            # remove all peaks that do not have an
            # average base pair q-value <=10^(-${q})
            awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} \
                $9>=q && $1!="ChrC" && $1!="ChrM"{print}' \
                ${macs2_out}/${chip_name}_R${t}_peaks.narrowPeak > \
                ${macs2_out}/noMask_qval${q}/${chip_name}_R${t}_peaks.narrowPeak

            # remove all peaks that overlap greenscreen
            bedtools intersect -v -wa \
                -a ${macs2_out}/noMask_qval${q}/${chip_name}_R${t}_peaks.narrowPeak \
                -b ${gs_regions} > \
                ${final_peaks}
        fi
    done


    # get pooled chip bams
    final_peaks="${macs2_out}/gsMask_qval${q}/${chip_name}_peaks.narrowPeak"
    if [[ ! -f $final_peaks ]]; then
        t_down_rep="${downsamp_bam_dir}/${chip_name}_R1.${bam_suffix}"
        t_pool_param="${t_down_rep}"
        for ((t=2; t<=${chip_nreps}; t++ )); do
            t_down_rep="${downsamp_bam_dir}/${chip_name}_R${t}.${bam_suffix}"
            t_pool_param="${t_pool_param} ${t_down_rep}"
        done
        echo "-t ${t_pool_param}"

        # call peaks on pooled ${chip_name}

        # run MACS2
        macs2 callpeak -t ${t_pool_param} \
            -c ${c_param} \
            -f BAM --keep-dup auto \
            --nomodel --extsize ${chip_fsize} -g 101274395 \
            --outdir ${macs2_out} -n ${chip_name}

        # remove all peaks that do not have an
        # average base pair q-value <=10^(-${q})
        awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} \
            $9>=q && $1!="ChrC" && $1!="ChrM"{print}' \
            ${macs2_out}/${chip_name}_peaks.narrowPeak > \
            ${macs2_out}/noMask_qval${q}/${chip_name}_peaks.narrowPeak

        # remove all peaks that overlap greenscreen
        bedtools intersect -v -wa \
            -a ${macs2_out}/noMask_qval${q}/${chip_name}_peaks.narrowPeak \
            -b ${gs_regions} > \
            ${macs2_out}/gsMask_qval${q}/${chip_name}_peaks.narrowPeak
    fi
done < meta/chip_controls_fragsize_nreps.csv
