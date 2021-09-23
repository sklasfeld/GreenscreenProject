#!/bin/bash

# output directory
out_dir="data/annotations/LFY_Jin_2021"

# annotate peaks upstream  to genes 
# with distance <=3000
ann_dist=-3000

# macs2 average basepair q-value threshold (log10)
q=10

# get lfy peak and summit bed file
lfy_peaks="data/macs2_out/chipPeaks/gsMask_qval${qthresh}/LFY_W_2021_peaks.broadPeak"
lfy_summits="data/macs2_out/chipPeaks/gsMask_qval${qthresh}/LFY_W_2021_summits.broadPeak"
awk -F"\t" 'BEGIN{OFS="\t"} \
	{print $1,$2+$10,$2+$10+1,$4,$5,$6,$7,$8,$9,"0"}' \
	${lfy_peaks} > ${lfy_summits}

# get gene bed files
gene_bed="meta/ArabidopsisGenome/Araport11_GFF3_nonHypothetical_proteinCoding_miRNA_genes.201606.bed"
lfy_gene_bed="meta/lfy_rna_diffExp/lfy_dependent_genes.bed"


nfield_gene_bed=`awk 'NR==1{print NF}' ${gene_bed}`
nfield_lfy_summits_bed=`awk 'NR==1{print NF}' ${lfy_summits}`

# annotate peaks to genes (round 1) 
#   * ignore peaks that are downstream of genes
#	* maximum upstream distance = 3000bp

bedtools closest -D a -id \
	-a ${gene_bed} \
	-b ${lfy_summits} | \
	awk -v minDist=${ann_dist} \
	-v nfield1=${nfield_gene_bed} \
	-v nfield2=${nfield_lfy_summits_bed} \
	-F"\t" \
	'BEGIN{OFS="\t"} \
	$NF>=minDist && $NF<=0{ \
		printf $4; \
		for(i=nfield1+1;i<=nfield1+nfield2;i++){printf "\t"$i} \
		printf "\n"}' \
	> ${out_dir}/lfy_summits_ann_round1.tsv


# get peaks not annotated in round 1
cut -f2- ${out_dir}/ann_round1.tsv \
	> ${out_dir}/lfy_summits_ann_round1.bed
comm -1 -3 \
	${out_dir}/lfy_summits_ann_round1.bed \
	${lfy_summits} \
	>  ${out_dir}/lfy_summits_notAnn_round1.bed

# annotate orphan peaks to genes (round 2)
#	* annotate to all lfy dependent genes within 10kb
round2_max_dist=10000

# get region +/-10kb of lfy summits
lfy_orphan_regions="${out_dir}/lfy_summitsPM${round2_max_dist}_notAnn_round1.bed"
awk -F"\t" -v max_dist=${round2_max_dist} \
	'BEGIN{OFS="\t"} \
	{if($2-max_dist) < 0{$2=0}else{$2=$2-max_dist} \
	print $1,$2,$3+max_dist-1,$4,$5,$6,$7,$8,$9,$10}' \
	${out_dir}/lfy_summits_notAnn_round1.bed \
	> {out_dir}/${lfy_orphan_regions}

bedtools intersect -wao \
	-a {out_dir}/${lfy_orphan_regions} \
	-b ${lfy_gene_bed} | \
	awk -v nfield1=${nfield_gene_bed} \
	-v nfield2=${nfield_lfy_summits_bed} \
	-F"\t" \
	'BEGIN{OFS="\t"} \
	{printf $(nfield2+4); \
	for(i=1;i<=nfield2;i++){printf "\t"$i} \
	printf "\n"} \
	> ${out_dir}/lfy_summits_ann_round2.tsv

cat \
	{out_dir}/lfy_summits_ann_round1.tsv \
	{out_dir}/lfy_summits_ann_round2.tsv | \
	awk -F"\t" \
	'BEGIN{OFS="\t"; \
	printf "gene_id", "chrom", "chipSummit_start";
	printf "\tchipSummit_stop", "chipPeak_name", "chipPeak_score"; \
	printf "\tchipPeak_strand", "chipPeak_signalValue"; \
	printf "\tchipPeak_pValue", "chipPeak_qValue"; \
	print "\tchipPeak_summitDist"} \
	{print}' > \
	{out_dir}/lfy_summits_ann.tsv
