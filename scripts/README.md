See supplemental methods for more detailed information.

# 1.3 Green screen Generation
1. scripts/import_raw_fasta_inputs.sh
	a. Purpose: batch imports raw FASTQ files
	b. Output Location: Current Working Directory
2. scripts/organize_raw_fasta_input.sh
	a. Purpose: batch organizes raw FASTQ files
	b. Output Location: fastq/raw
3. scripts/fastqc_raw_inputs.sh
	a. Purpose: batch FASTQC on raw FASTQ files
	b. Output Location: data/FASTQC/raw
4. scripts/trimmomatic_inputs.sh 
	a. Purpose:
		i. batch TRIMMOMATIC on raw FASTQ files
		ii. batch FASTQC on TRIMMOMATIC FASTQ output
	b. Output Location: 
		i. fastq files: fastq/trimmed
		ii. fastqc files: data/FASTQC/trimmed 
5. scripts/mapped_inputs.sh
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
6. scripts/predictd_inputs.sh
	a. Purpose: batch measure strand bias with macs2 predictd
	b. Steps:
		i. Run macs2 predictd
		ii. Run Rscript output from macs2 predictd
		iii. Put the pdf output from the Rscript into the correct directory
	c. Output Location: data/plotStrandBias/inputs
7. scripts/macs2_callpeaks_inputs.sh
	a. Purpose: batch call peaks in each input with MACS2
	b. Meta Requirement: meta/input_readsizes.tsv
	c. Steps
		i. Run macs2 callpeaks
		ii. Filter peaks by average log10 q-value (10)
	d. Output Location: data/macs2_out/inputControls
8. scripts/generate_20input_greenscreenBed.sh
	a. Purpose: batch generate greenscreen from peaks called in 20 individual inputs
	b. Steps:
		i. Concatenate peaks
		ii. Merge all overlapping regions and regions within (1000) bp apart
		iii. Filter out regions that are called in less than (10) distinct samples
	c. Output Location: data/arabidopsis_greenscreen_20inputs.bed
# 1.4 Calling ChIP-Seq Peaks
1. scripts/import_raw_fasta_publishedChIPsAndControls.sh 
	a. Purpose: batch imports raw FASTQ files
	b. Output Location: Current Working Directory
2. scripts/organize_raw_fasta_publishedChIPsAndControls.sh
	a. Purpose: batch organizes raw FASTQ files
	b. Output Location: fastq/raw
3. scripts/fastqc_raw_fasta_publishedChIPsAndControls.sh
	a. Purpose: batch FASTQC on raw FASTQ files
	b. Meta Requirement: meta/chip_controls_sampleIDs.list
	c. Output Location: data/FASTQC/raw
4. scripts/trimmomatic_publishedChIPsAndControls.sh
	a. Purpose:
		 i. batch TRIMMOMATIC on raw FASTQ files
		 ii. batch FASTQC on TRIMMOMATIC FASTQ output
	b. Output Location: 
		i. fastq files: fastq/trimmed
		ii. fastqc files: data/FASTQC/trimmed 
5. scripts/mapped_publishedChIPsAndControls.sh
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
6. scripts/predictd_publishedChIPsAndControls.sh
	a. Purpose: batch measure strand bias with macs2 predictd
	b. Meta Requirement: meta/chip_controls_sampleIDs.list
	c. Steps:
		i. Run macs2 predictd
		ii. Run Rscript output from macs2 predictd
		iii. Put the pdf output from the Rscript into the correct directory
	d. Output Location: data/plotStrandBias/chip
7. ChIPQC.R
	a. Purpose: returns ChIP-Seq quality control metrics
8. scripts/chipReplicates_BamToBigwig.sh
	a. Purpose: generate IGV bigwig files for ChIP-Seq replicates
	b. Meta Requirement: meta/chipFragSize.csv
	c. Steps
		i. Convert BAM-format to BED-format
		ii. Sort the BED file
		iii. Extend the reads in the BED file
		iv. Normalize the signal and output BEDGRAPH-format
		v. compress BEDGRAPH to BIGWIG format
	d. Output Location: mapped/chip
9. scripts/chipPooledBigwig.sh
	a. Purpose: generate IGV bigwig files for ChIP-Seq samples (average of the replicates)
	b. Steps
		i. Get matrix of signal for each replicate
		ii. Find average signals
		iii. Report the average signal in BIGWIG format
	c. Output Location: mapped/chip
10. scripts/downSampleBam_publishedChIPsAndControls.sh
	a. Purpose: downsample each replicate to match the sampling depth of the replicate with the lowest sampling depth.
	b. Output Location: mapped/chip/downsample
11. scripts/macs2_callpeaks_publishedChIPs.sh
	a. Purpose: call peaks on published ChIPs
	b. Steps
		i. Run macs2 callpeaks
		ii. Filter peaks by average log10 q-value (10)
	c. Output Location: data/macs2_out/chipPeaks/qval10
