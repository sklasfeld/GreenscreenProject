#!/bin/bash

raw_bam_dir="mapped/chip"
downsamp_bam_dir="${raw_bam_dir}/downsample"
bam_suffix="dupmark.sorted.bam"
macs2_out="data/macs2_out/chipPeaks"
# average basepair q-value threshold (log10)
q=10

# make macs2 output directory
mkdir -p ${macs_out}/qval${q}(

callpeak_cmd (){
        treatment="$1"
        control="$2"
        esize="$3"
        odir="$4"
        oname="$5"
        qthresh="$6"

        # run MACS2
        macs2 callpeak -t ${treatment} -c ${control} -f BAM --keep-dup auto \
                --nomodel --extsize ${esize} -g 101274395 \
                --outdir ${odir} -n ${oname}

        # remove all peaks that do not have an
        # average base pair q-value <=10^(-${q})
        awk -F"\t" -v q=${qthresh} 'BEGIN{OFS="\t"} \
                $9>=q && $1!="ChrC" && $1!="ChrM"{print}' \
                ${odir}/${oname}_peaks.broadPeak > \
                ${odir}/qval${qthresh}/${oname}_peaks.broadPeak
}

while read line; do
        chip_name=$(echo $line | cut -d "," -f1)
        chip_nreps=$(echo $line | cut -d "," -f2)
        chip_fsize=$(echo $line | cut -d "," -f2)
        cntl_name=$(echo $line | cut -d "," -f4)
        cntl_nreps=$(echo $line | cut -d "," -f5)
        
        c_param="${downsamp_bam_dir}/${cntl_name}_R1.${bam_suffix}"
        for ((c=2; c<=${cntl_nreps}; c++ )); do
                c_param="${c_param} ${downsamp_bam_dir}/${cntl_name}_R${c}.${bam_suffix}"
        done
        echo "-c ${c_param}"

        # call peaks on ${chip_name}_R1
        t_raw_rep="${raw_bam_dir}/${chip_name}_R1.${bam_suffix}"
        callpeak_cmd ${t_raw_rep} ${c_param} ${chip_fsize} \
                ${macs2_out} ${chip_name}_R1 ${q}

        # get pooled chip bams
        t_down_rep="${downsamp_bam_dir}/${chip_name}_R1.${bam_suffix}"
        t_pool_param="${t_down_rep}"
        for ((t=2; t<=${chip_nreps}; t++ )); do
                # call peaks on ${chip_name}_R${t}
                t_raw_rep="${downsamp_bam_dir}/${chip_name}_R${t}.${bam_suffix}"
                callpeak_cmd ${t_raw_rep} ${c_param} ${chip_fsize} \
                        ${macs2_out} ${chip_name}_R${t} ${q}

                # get pooled chip bams
                t_pool_param="${t_pool_param} ${t_raw_rep}"
        done
        echo "-t ${t_pool_param}"
        
        # call peaks on pooled ${chip_name}
        callpeak_cmd ${t_pool_param} ${c_param} ${chip_fsize} \
                ${macs2_out} ${chip_name} ${q}

done < meta/chip_controls_nreps.csv

