#!/bin/bash

# directory to output MACS2 results
macs_out="data/macs2_out/inputControls"
# average basepair q-value threshold (log10)
q=10

# make macs2 output directory
mkdir -p ${macs_out}/qval${q}

while read line; do
    x=$(echo $line | cut -d "," -f1)
    readsize=$(echo $line | cut -d "," -f2)

    # call peaks with MACS2
    macs2 callpeak \
        -t mapped/input/input${x}.dupmark.sorted.bam \
        -f BAM --keep-dup auto --nomodel \
        --extsize ${readsize} --broad --nolambda \
        -g 101274395 -n input${x} \
        --outdir ${macs_out}

    # remove all peaks that do not have an
    # average base pair q-value <=10^(-${q})
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} \
        $9>=q && $1!="ChrC" && $1!="ChrM"{print}' \
        ${macs_out}/input${x}_peaks.broadPeak > \
        ${macs_out}/qval${q}/input${x}_peaks.broadPeak
done < meta/input_readsizes.csv
