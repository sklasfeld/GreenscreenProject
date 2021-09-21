#!/bin/bash

greenscreen_regions="data/macs2_out/inputControls/qval10/gs_merge5000bp_call10_20inputs.txt"

# 1. mask reads that overlap Greenscreen

# create output directory for reads
# that do and do not overlap
# greenscreen regions
mkdir -p mapped/chip/greenscreen_mask
mkdir -p mapped/chip/greenscreen_regions


while read line; do
    samp=`echo $line | cut -d "," -f1`
    if [[ "$samp" != "SampleID" ]]; then
	# split reads that do and do not overlap
	# greenscreen
    	index_bam="mapped/chip/greenscreen_mask/${samp}.dupmark.sorted.bam.bai"
    	if [[ ! -f "$index_bam" ]]; then
    	    samtools view -b -L ${greenscreen_regions} \
                -o mapped/chip/greenscreen_regions/${samp}.dupmark.sorted.bam \
                -U mapped/chip/greenscreen_mask/${samp}.dupmark.sorted.bam \
                mapped/chip/${samp}.dupmark.sorted.bam
    	    samtools index mapped/chip/greenscreen_mask/${samp}.dupmark.sorted.bam
        fi
    fi
done < meta/gsMaskReads_publishedChIPsAndControls_sampleSheet.csv

# 2. Run ChIPQC

# create an output directory for ChIPQC results
mkdir -p data/ChIPQCreport/chip_gsMask_wiDups

Rscript scripts/ChIPQC.R \
    -f Replicate -c Factor \
    --indivReports -g Araport11 \
    -c Chr1 Chr2 Chr3 Chr4 Chr5 \
    -a meta/ArabidopsisGenome/Araport11_GFF3_genes_transposons.201606.gff \
    -s meta/ArabidopsisGenome/TAIR10_chr_count.txt \
    gsMaskReads_publishedChIPsAndControls_sampleSheet.csv \
    data/ChIPQCreport/chip_gsMask_wiDups
