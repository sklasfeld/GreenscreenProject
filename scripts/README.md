See supplemental methods for more detailed information.

# 1.3 Green screen Generation
1. `scripts/import_raw_fasta_inputs.sh`
	
	a. Purpose: batch imports raw FASTQ files
	
	b. Output Location: Current Working Directory

2. `scripts/organize_raw_fasta_input.sh`
	
	a. Purpose: batch organizes raw FASTQ files
	
	b. Output Location: fastq/raw

3. `scripts/fastqc_raw_inputs.sh`
	
	a. Purpose: batch FASTQC on raw FASTQ files
	
	b. Output Location: data/FASTQC/raw

4. `scripts/trimmomatic_inputs.sh`
	
	a. Purpose:
		
		i. batch TRIMMOMATIC on raw FASTQ files
		
		ii. batch FASTQC on TRIMMOMATIC FASTQ output
	
	b. Output Location: 
		
		i. fastq files: fastq/trimmed
		
		ii. fastqc files: data/FASTQC/trimmed 

5. `scripts/mapped_inputs.sh`

	a. Purpose: batch generate, sort, and index mapped reads with MAPQ >=30
	
	b. Steps:
	
		i. Run bowtie2 on trimmed FASTQ files
		
		ii. Sort the reads
		
		iii. Index the bam file
		
		iv. Remove reads without MAPQ>=30
		
		v. Index filtered reads
		
		vi. Mark duplicates with picard
		
		vii. Sort reads after marking the duplicates
		
		viii. Index the sorted reads
		
	c. Output Location: mapped/input

7. `scripts/translateUniCodeFile.py`

	a. Purpose: translate files that contain percent unicode symbols

8. `scripts/measureContigLengthFromFasta.sh`

	a. Purpose: Calculate genome contig sizes

9. `scripts/ChIPQC.R`
	
	a. Purpose: returns ChIP-Seq quality control metrics

10. `scripts/input_BamToBigwig.sh`

	a. Purpose: convert mapped reads into compressed bigwig files 
	to use for IGV visualization of the signals

	b. Output Location: data/

11. `scripts/macs2_callpeaks_inputs.sh`
	
	a. Purpose: batch call peaks in each input with MACS2
	
	b. Meta Requirement: meta/input_readsizes.tsv
	
	c. Steps
		
		i. Run macs2 callpeaks
		
		ii. Filter peaks by average log10 q-value (10)
	
	d. Output Location: data/macs2_out/inputControls
	
12. `scripts/generate_20input_greenscreenBed.sh`
	
	a. Purpose: batch generate greenscreen from peaks called in 20 individual inputs
	
	b. Steps:
	
		i. Concatenate peaks
		
		ii. Merge all overlapping regions and regions within (1000) bp apart
		
		iii. Filter out regions that are called in less than (10) distinct samples
	
	c. Output Location: data/arabidopsis_greenscreen_20inputs.bed

# 1.4 Calling ChIP-Seq Peaks

1. `scripts/import_raw_fasta_publishedChIPsAndControls.sh`
	
	a. Purpose: batch imports raw FASTQ files
	
	b. Output Location: Current Working Directory

2. `scripts/organize_raw_fasta_publishedChIPsAndControls.sh`
	
	a. Purpose: batch organizes raw FASTQ files
	
	b. Output Location: fastq/raw

3. `scripts/fastqc_raw_fasta_publishedChIPsAndControls.sh`
	
	a. Purpose: batch FASTQC on raw FASTQ files
	
	b. Meta Requirement: meta/chip_controls_sampleIDs.list
	
	c. Output Location: data/FASTQC/raw

4. `scripts/trimmomatic_publishedChIPsAndControls.sh`
	
	a. Purpose:
		 
		 i. batch TRIMMOMATIC on raw FASTQ files
		 
		 ii. batch FASTQC on TRIMMOMATIC FASTQ output
	
	b. Output Location: 
		
		i. fastq files: fastq/trimmed
		
		ii. fastqc files: data/FASTQC/trimmed 
		
5. `scripts/mapped_publishedChIPsAndControls.sh`
	
	a. Purpose: batch generate, sort, and index mapped reads with MAPQ >=30
	
	b. Meta Requirement: meta/chip_controls_sampleIDs.list
	
	c. Steps:
		
		i. Run bowtie2 on trimmed FASTQ files
		
		ii. Sort the reads
		
		iii. Index the bam file
		
		iv. Remove reads without MAPQ>=30
		
		v. Index filtered reads
		
		vi. Mark duplicates with picard
		
		vii. Sort reads after marking the duplicates
		
		viii. Index the sorted reads
	
	d. Output Location: mapped/chip

6. `scripts/ChIPQC.R` (repeat)
	
	a. Purpose: returns ChIP-Seq quality control metrics

7. `scripts/ChIPQC_gsMaskReads_publishedChIPsAndControls.sh`
	
	a. Purpose: call ChIPQC library on each ChIP-Seq experiment
	
	b. Meta Requirement: meta/gsMaskReads_publishedChIPsAndControls_sampleSheet.csv
	
	c. Steps:
		
		i. mask reads that overlap greenscreen in each ChIP-Seq and control experiment
		
		ii. Run ChIPQC library commands on each sample
		
	d. Output Location: data/ChIPQCreport/chip_gsMask_wiDups

8. `scripts/chipReplicates_BamToBigwig.sh`
	
	a. Purpose: generate IGV bigwig files for ChIP-Seq replicates
	
	b. Meta Requirement: meta/chipFragSize.csv
	
	c. Steps
		
		i. Convert BAM-format to BED-format
		
		ii. Sort the BED file
		
		iii. Extend the reads in the BED file
		
		iv. Normalize the signal and output BEDGRAPH-format
		
		v. compress BEDGRAPH to BIGWIG format
	
	d. Output Location: mapped/chip

9. `scripts/chipPooledBigwig.sh`
	
	a. Purpose: generate IGV bigwig files for ChIP-Seq samples (average of the replicates)
	
	b. Steps
		
		i. Get matrix of signal for each replicate
		
		ii. Find average signals
		
		iii. Report the average signal in BIGWIG format
	
	c. Output Location: mapped/chip

10. `scripts/downSampleBam_publishedChIPsAndControls.sh`
	
	a. Purpose: downsample each replicate to match the sampling depth of the replicate with the lowest sampling depth.
	
	b. Output Location: mapped/chip/downsample

11. `scripts/macs2_callpeaks_publishedChIPs.sh`
	
	a. Purpose: call peaks on published ChIPs
	
	b. Steps
		
		i. Run macs2 callpeaks
		
		ii. Filter peaks by average log10 q-value (10)

		iii. Remove peaks that overlap greenscreen regions
	
	c. Output Location: data/macs2_out/chipPeaks/gsMask_qval10

# 1.5 ChIP-Seq Analysis
## Annotate LFY ChIP-Seq Peaks

1. `scripts/generateGeneBedFile.sh`

	a. Purpose: Given the Arabidopsis genome annotation, export ALL protein-coding
	and miRNA gene locations in BED-format

	b. Output Location: meta/ArabidopsisGenome/Araport11_GFF3_nonHypothetical_proteinCoding_miRNA_genes.201606.bed

2.`scripts/generateLFYDependentGeneBedFile.sh`

	a. Purpose: Given the Arabidopsis genome annotation, export protein-coding
	and miRNA gene locations in BED-format of genes that were found differentially
	expressed after LFY gene expression

	b. Output Location: scripts/generateLFYDependentGeneBedFile.sh

3. `scripts/annotate_LFY_peaks.sh`

	a. Purpose: Annotate LFY peaks from Jin 2021

	b. Steps

		i. First annotate all peaks to the closest genes of which peaks are upstream of genes. Peaks should be at most 3kb upstream of a genes. (round1 annotation)

		ii. Of peaks not annotated in the previous step (orphan peaks), annotate to genes within 10kb that are differentially expressed in the presence of LFY (round2 annotation)
	
	c. Output Location: data/annotations/LFY_Jin_2021/lfy_summits_ann.tsv
	
## Compare ChIP-Seq samples from different publications

1. `scripts/merge_chip_peaks.sh`

	a. Purpose: peaks called in each ChIP-seq experiment replicate are concatonated and merged

	b. Output Location: data/macs2_out/chipPeaks/gsMask_qval10/ChIPseq_Peaks.merged.bed

2. `scripts/coverage_bed_matrix.py`

	a. Purpose: Given specific regions and multiple bigwig files, 
	generate a signal matrix where each row represents a specific
	region and each column represents a sequencing experiment

	b. Example command: 


``` {.bash frame="lines" breaklines=""}
python3 scripts/coverage_bed_matrix.py \
    meta/chip_trueRep_bigwigs.csv \
    data/macs2_out/chipPeaks/gsMask_qval10/ChIPseq_Peaks.merged.bed \
    -o data/plotCorrelation \
    -m coverage_matrix_trueRep_peaks_merged.csv
```

3. `scripts/readCorrelationPlot.py`

	a. Purpose: Given a signal matrix where each row represents a specific region and each column represents a sequencing experiment, calculate
	a correlation coefficient between samples and display these values in a heatmap sorted using hierarchical clustering. Optional: calculate 
	rand-index values given expected clusters

	b. Example command:

``` {.bash frame="lines" breaklines=""}
python3 scripts/readCorrelationPlot.py \
    data/plotCorrelation/coverage_matrix_trueRep_peaks_merged.csv \
    data/plotCorrelation/trueRep_peaks_merged_heatmap.png \
    -lm ward --plot_numbers -k 2 -ri \
    -sl meta/chip_trueReps_colorshapeLabels.csv \
    -cf meta/chip_trueReps_expectedCluster.csv
```