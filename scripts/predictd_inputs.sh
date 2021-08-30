#!/bin/bash

# create output directories
plot_out="data/plotStrandBias/inputs"
script_out=${plot_out}/scripts

input_list=("A" "B" "C" "D" "E"
        "F" "G" "H" "I" "J" "K"
        "L" "M" "N" "O" "P" "Q"
        "R" "S" "T")

for x in "${input_list[@]}"; do
    # run macs2 predictd
    macs2 predictd -g 101274395 \
        --outdir ${script_out} \
        --rfile predictd_input${x}.R \
        -i mapped/input/input${x}.dupmark.sorted.bam
    # run Rscript output from macs2 predict d
    Rscript ${script_out}/predictd_input${x}.R
    # put the pdf output from the Rscript into the correct directory
    mv predictd_input${x}.R_model.pdf ${plot_out}
done
