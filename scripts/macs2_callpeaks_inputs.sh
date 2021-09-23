#!/bin/bash

# directory to output MACS2 results
macs_out="data/macs2_out/inputControls"

# make macs2 output directory
mkdir -p ${macs_out}

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

done < meta/input_readsizes.csv
