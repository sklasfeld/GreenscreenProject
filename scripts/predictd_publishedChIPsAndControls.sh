#!/bin/bash

greenscreen_regions="data/macs2_out/inputControls/qval10/concat_input_peaks.broadPeak"

# create output directories
plot_out="data/plotStrandBias/chip"
script_out=${plot_out}/scripts

# create directory for reads
# that do and do not overlap
# greenscreen regions
mkdir -p mapped/chip/greenscreen_mask
mkdir -p mapped/chip/greenscreen_regions

while read samp; do
    # analyze all reads (no Mask)
    #   run macs2 predictd
    macs2 predictd -g 101274395 \
        --outdir ${script_out} \
        --rfile predictd_noMask_${samp}.R \
        -i mapped/chip/${samp}.dupmark.sorted.bam
    #   run Rscript output from macs2 predict d
    Rscript ${script_out}/predictd_noMask_${samp}.R
    #   put the pdf output from the Rscript into the correct directory
    mv predictd_noMask_${samp}.R_model.pdf ${plot_out}

    # split reads that do and do not overlap 
    # greenscreen
    samtools view -b -L ${greenscreen_regions} \
        -o mapped/chip/greenscreen_regions/${samp}.dupmark.sorted.bam \
        -U mapped/chip/greenscreen_mask/${samp}.dupmark.sorted.bam \
        mapped/chip/${samp}.dupmark.sorted.bam
    samtools index mapped/chip/greenscreen_mask/${samp}.dupmark.sorted.bam

    # analyze reads that do not overlay
    # green screen regions
    #   run macs2 predictd
    macs2 predictd -g 101274395 \
        --outdir ${script_out} \
        --rfile predictd_gsMask_${samp}.R \
        -i mapped/chip/greenscreen_mask/${samp}.dupmark.sorted.bam
    #   run Rscript output from macs2 predict d
    Rscript ${script_out}/predictd_gsMask_${samp}.R
    #   put the pdf output from the Rscript into the correct directory
    mv predictd_gsMask_${samp}.R_model.pdf ${plot_out}
done < meta/chip_controls_sampleIDs.list
