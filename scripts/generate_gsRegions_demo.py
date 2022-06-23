#!/usr/bin/env python

import os
import sys
import subprocess
import tempfile
import pandas as pd
import numpy as np

# function to add a directory
# only if the directory does not
# exist
def add_dir(dir_path):
	if not os.path.isdir(dir_path):
		os.makedirs(dir_path)

# function to measure number of reads
# in sam/bam file
def bamsize(bamfile) :
	cmd=("samtools view -c %s" % bamfile)
	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
	return(int(proc.communicate()[0].decode("utf-8").strip()))

def isempty(file_path):
	if os.path.isfile(file_path):
		if (file_path.split(".")[-1] == "bam"
			or file_path.split(".")[-1] == "sam"):
				if bamsize(file_path) == 0:
					return(True)
				else:
					return(False)
		else:
			if os.stat(file_path).st_size == 0:
				return (True)
			else:
				return(False)
	else:
		return(True)

# function to run and print BASH commands
def cmd(foo, verbose=False, output_f=""):
	run_cmd=True
	if len(output_f) > 0:
		if not isempty(output_f):
			run_cmd=False
	if run_cmd:
		if verbose:
			sys.stdout.write("> CMD: %s\n" % foo)
			sys.stdout.flush()
		os.system(foo)

# Hard-Coded Variables

# * meta input information
input_list=["inputA_demo", "inputB_demo", "inputC_demo"]
readsize_dic={"inputA_demo":100,
	"inputB_demo":80,
	"inputC_demo":76}

# * read filtering parameters
samtools_flag_filter = 772
samtools_mapq = 30

# * normalization settings for bigwig format
normalizeby = 10000000 # scaling factor

# * greenscreen settings
qval=10
merge_distance=5000
distinct_ninputs=2

# * other settings
verbose = True
chromosomes=["Chr2"]
genome_prefix=("TAIR10_%s" % ("_".join(chromosomes)))

# * paths
meta_dir = "demo/meta"
data_dir = "demo/data"
genome_meta_dir = ("%s/ArabidopsisGenome" % meta_dir)
genome_fasta = ("%s/%s.fasta" % (genome_meta_dir, genome_prefix))
gff_file = ("%s/Araport11_GFF3_genes_transposons.UPDATED.201606.Chr2.gff"
	% genome_meta_dir)
fastq_raw_dir = ("%s/fastq/raw" % data_dir)
fastqc_raw_dir = ("%s/FASTQC/trimmed" % data_dir)
trimmomatic_install_dir = "/usr/src/Trimmomatic-0.39"
fastq_trim_dir = ("%s/fastq/trimmed" % data_dir)
fastqc_trim_dir = ("%s/FASTQC/trimmed" % data_dir)
bowtie2_index_dir = ("%s/bowtie2_genome_dir" % genome_meta_dir)
bowtie2_index_prefix  =  ("%s/%s" % (bowtie2_index_dir, genome_prefix))
chipqc_dir = ("%s/ChIPQCreport/inputs_noMask" % data_dir)
picard_path = "/usr/src/picard/build/libs"
mapped_dir = ("%s/mapped/input" % data_dir)
macs_dir = ("%s/macs2_out/inputControls" % data_dir)
qval_dir = ("%s/qval%i" % (macs_dir, qval))
concat_file = ("%s/concat_input_peaks.broadPeak" % qval_dir)
merge_file = ("%s/merge%ibp_demoInputs.txt" % (qval_dir,merge_distance))
final_greenscreen = ("%s/gs_merge%ibp_call%i_demoInputs.broadPeak" 
	% (qval_dir, merge_distance,distinct_ninputs))



# 1. count genome chromosome sizes
chrom_sizes = (
	"%s/%s_count.txt"
	% (genome_meta_dir, genome_prefix))
chrom_sizes_cmd = (
	("bash scripts/measureContigLengthFromFasta.sh %s %s")
	% (genome_fasta,chrom_sizes))
cmd(chrom_sizes_cmd, verbose, chrom_sizes)
# * check output
if isempty(chrom_sizes):
	sys.exit("Error: Custom script " +
		"`scripts/measureContigLengthFromFasta.sh` " +
		"did not generate table as expected\n")

# 2. build genome bowtie2 index
add_dir(bowtie2_index_dir)
if not os.path.isdir(bowtie2_index_dir):
	bowtie2_index_cmd = (
		"bowtie2-build %s %s" %
		(genome_fasta, bowtie2_index_prefix))
	cmd(bowtie2_index_cmd, verbose)
elif len(os.listdir(bowtie2_index_dir)) == 0:
	bowtie2_index_cmd = (
		"bowtie2-build %s %s" %
		(genome_fasta, bowtie2_index_prefix))
	cmd(bowtie2_index_cmd, verbose)
bowtie2_index_err=("Error: Bowtie2 did not build genome index\n")
# * check output
if not os.path.isdir(bowtie2_index_dir):
	sys.exit(bowtie2_index_err)
if len(os.listdir(bowtie2_index_dir)) == 0:
	sys.exit(bowtie2_index_err)


round=1
ninputs=len(input_list)
broadPeak_l=[]
for samp in input_list:
	raw_fastq = ("%s/%s.fq.gz" % (fastq_raw_dir,samp))
	
	# 3. run fastqc on raw fastq
	add_dir(fastqc_raw_dir)
	fastqc_raw_output = ("%s/%s_fastqc.html" %
		(fastqc_raw_dir, samp))
	fastqc_raw_cmd = (("fastqc -o %s %s") % 
		(fastqc_raw_dir, raw_fastq))
	cmd(fastqc_raw_cmd, verbose, fastqc_raw_output)
	# * check output
	
	if isempty(fastqc_raw_output):
		sys.exit("Error: FASTQC output for %s reads were not generated\n"
			% samp)
	
	# 4. run trimming
	add_dir("adapters")
	cp_adapters_cmd = ("cp %s/adapters/TruSeq3-SE.fa adapters"
		% (trimmomatic_install_dir))
	cmd(cp_adapters_cmd, verbose, "adapters/TruSeq3-SE.fa")
	
	add_dir(fastq_trim_dir)
	trim_fastq = ("%s/%s.trimmed.fq.gz" % (fastq_trim_dir,samp))
	trim_cmd = (("java -jar %s/trimmomatic-0.39.jar "+
		"SE -phred33 %s %s " + 
		"ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 " + 
		"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33") % 
		(trimmomatic_install_dir, raw_fastq, trim_fastq))
	cmd(trim_cmd, verbose, trim_fastq)
	# * check output
	if isempty(trim_fastq):
		sys.exit("Error: Trimmomatic output for %s was not generated\n"
			% samp)

	# 5. run fastqc on trimmed fastq
	add_dir(fastqc_trim_dir)
	fastqc_trim_output = ("%s/%s.trimmed_fastqc.html" %
		(fastqc_trim_dir, samp))
	fastqc_trim_cmd = ("fastqc -o %s %s"
		%(fastqc_trim_dir, trim_fastq))
	cmd(fastqc_trim_cmd, verbose, fastqc_trim_output)
	# * check output
	if isempty(fastqc_trim_output):
		sys.exit("Error: FASTQC output for %s trimmed reads were not generated\n"
			% samp)
	
	# 6. map reads to genome with bowtie2
	add_dir(mapped_dir)
	bowtie2_mapped_f = ("%s/%s.sam" %
		(mapped_dir, samp))
	bowtie2_map_cmd = (("bowtie2  --phred33 -q -x %s -S %s %s") % 
		(bowtie2_index_prefix, bowtie2_mapped_f, trim_fastq))
	cmd(bowtie2_map_cmd, verbose, bowtie2_mapped_f)
	# * check output
	bowtie2_mapped_file_err = (("Error: bowtie2 mapping output for %s was " +
		"not generated\n") % samp)
	if isempty(bowtie2_mapped_f):
		sys.exit(bowtie2_mapped_file_err)


	# 7. sort the reads
	bowtie2_mapped_sorted_f = ("%s/%s.bam" % (mapped_dir, samp))
	bowtie2_mapped_sorted_cmd = (
		"samtools sort -o %s %s" % 
		(bowtie2_mapped_sorted_f, bowtie2_mapped_f))
	cmd(bowtie2_mapped_sorted_cmd, verbose, bowtie2_mapped_sorted_f)
	# * check output
	bowtie2_mapped_sorted_file_err = (("Error: bowtie2 mapping output for " +
		"%s was not sorted\n") % samp)
	if isempty(bowtie2_mapped_sorted_f):
		sys.exit(bowtie2_mapped_sorted_file_err)

	# 8. index the bam file
	bowtie2_mapped_sorted_idx = ("%s.bai" % bowtie2_mapped_sorted_f)
	index_sorted_bowtie2_mapped_cmd = (
		"samtools index %s" % bowtie2_mapped_sorted_f)
	cmd(index_sorted_bowtie2_mapped_cmd, verbose, 
		bowtie2_mapped_sorted_idx)
	# * check output
	bowtie2_mapped_sorted_idx_err = (("Error: sorted bowtie2 mapping " +
		"output for %s was not indexed\n") % samp)
	if isempty(bowtie2_mapped_sorted_idx):
		sys.exit(bowtie2_mapped_sorted_idx_err)

	# 9. remove reads without MAPQ>=30 and without sam flags 772
	bowtie2_filtered_f = ("%s/%s.filter.bam" % (mapped_dir, samp))
	bowtie2_samtools_filter_cmd = (("samtools view -F %i " + 
		"-q %i -b %s -o %s %s") % (samtools_flag_filter, 
		samtools_mapq, bowtie2_mapped_sorted_f, bowtie2_filtered_f,
		"_".join(chromosomes)))
	cmd(bowtie2_samtools_filter_cmd, verbose, bowtie2_filtered_f)
	# * check output
	bowtie2_filter_file_err = (("Error: bowtie2 mapped output for " +
		"%s could not be filtered (MAPQ >= %i and only include " +
		"reads with none of the FLAGS in %i present)\n") % 
		(samp, samtools_mapq, samtools_flag_filter))
	if isempty(bowtie2_filtered_f):
		sys.exit(bowtie2_filter_file_err)

	# 10. sort filtered reads
	bowtie2_sorted_filtered_f = ("%s/%s.filter.sorted.bam" % (mapped_dir, samp))
	bowtie2_sort_filtered_cmd = ("samtools sort %s -o %s" 
		% (bowtie2_filtered_f, bowtie2_sorted_filtered_f))
	cmd(bowtie2_sort_filtered_cmd, verbose,
		bowtie2_sorted_filtered_f)
	# * check output
	bowtie2_sorted_filtered_file_err = (("Error: bowtie2 filtered output for " +
		"%s was not sorted\n") % samp)
	if isempty(bowtie2_sorted_filtered_f):
		sys.exit(bowtie2_sorted_filtered_file_err)

	# 11. generate index for filtered reads
	bowtie2_sorted_filtered_idx = ("%s.bai" % bowtie2_sorted_filtered_f)
	index_sorted_filtered_bowtie2_cmd = (
		"samtools index %s" % bowtie2_sorted_filtered_f)
	cmd(index_sorted_filtered_bowtie2_cmd, verbose, 
		bowtie2_sorted_filtered_idx)
	# * check output
	bowtie2_sorted_filtered_idx_err = (("Error: sorted filtered mapped reads " +
		"from bowtie2 for %s was not indexed\n") % samp)
	if isempty(bowtie2_sorted_filtered_idx):
		sys.exit(bowtie2_sorted_filtered_idx_err)

	# 12. mark duplicates with picard
	picard_out = ("%s/%s.dupmark.bam" % (mapped_dir, samp))
	markdups_cmd = (("java -jar %s/picard.jar MarkDuplicates " +
		"-I %s -O %s " +
		"-M %s/%s.dup.qc " +
		"-VALIDATION_STRINGENCY LENIENT " +
		"-REMOVE_DUPLICATES false " +
		"-ASSUME_SORTED true") %
		(picard_path, bowtie2_sorted_filtered_f, picard_out,
		mapped_dir, samp))
	cmd(markdups_cmd, verbose, picard_out)
	# * check output
	markdups_err = (("Error: duplicate were not marked for %s\n") % samp)
	if isempty(picard_out):
		sys.exit(markdups_err)

	# 13. sort reads after marking the duplicates
	markdups_sorted_f = ("%s/%s.dupmark.sorted.bam" % (mapped_dir, samp))
	markdups_sort_cmd = ("samtools sort %s -o %s" 
		% (picard_out, markdups_sorted_f))
	cmd(markdups_sort_cmd, verbose, markdups_sorted_f)
	# * check output
	markdups_sorted_file_err = (("Error: picard markdups output for " +
		"%s was not sorted\n") % samp)
	if isempty(markdups_sorted_f):
		sys.exit(markdups_sorted_file_err)


	# 14. Generate index for the sorted marked reads
	markdups_sorted_idx = ("%s.bai" % markdups_sorted_f)
	index_markdups_cmd = (
		"samtools index %s" % markdups_sorted_f)
	cmd(index_markdups_cmd, verbose, markdups_sorted_idx)
	# * check output
	markdups_idx_err = ("Error: picard output for %s was not indexed\n" % samp)
	if isempty(markdups_sorted_idx):
		sys.exit(markdups_idx_err)


	# 15. Generate sample sheet for ChIPQC code
	samplesheet=("%s/chipqc_sampleSheet_%s.csv" % (meta_dir, samp))
	samplesheet_handle = open(samplesheet, "w")
	samplesheet_handle.write("SampleID,bamReads\n")
	samplesheet_handle.write("%s,%s/%s.dupmark.sorted.bam\n" % 
		(samp, mapped_dir,samp))
	samplesheet_handle.close()

	# 15. Run ChIPQC
	add_dir(chipqc_dir)
	chipqc_out = (("%s/%s/ChIPQC.html") % (chipqc_dir, samp))
	chipqc_r_cmd = (("Rscript scripts/ChIPQC.R --indivReports -g Araport11 " +
		"-c %s -a %s -s %s %s %s") % 
		("_".join(chromosomes), gff_file, 
			chrom_sizes, samplesheet, chipqc_dir))
	cmd(chipqc_r_cmd, verbose, chipqc_out)

	# 16. normalize signal and output BEDGRAPH
	totreads= bamsize(markdups_sorted_f)
	scaling=normalizeby / totreads
	
	samp_bg = ("%s/%s.bg" % (mapped_dir, samp))
	if isempty(samp_bg):
		temp_bg = tempfile.NamedTemporaryFile()
		scale_bam = ("genomeCoverageBed -ibam %s -bg -scale %f > %s" 
			% (markdups_sorted_f, scaling, temp_bg.name))
		cmd(scale_bam, verbose)
	
		bg_handle = open(samp_bg, "w")
		with open(temp_bg.name) as temp_bg_handle:
			for line in temp_bg_handle:
				row = line.split("\t")
				row[3] = ("%.2f" % float(row[3]))
				bg_handle.write("\t".join(row))
		bg_handle.close()
		temp_bg.close()
		# * check output
		bedgraph_err = ("Error: %s bedgraph file was not generated.\n" % samp)
		if isempty(samp_bg):
			sys.exit(bedgraph_err)

	# 17. BEDGRAPH to BIGWIG FORMAT
	samp_bw = ("%s/%s.bw" % (mapped_dir, samp))
	if not os.path.isfile(samp_bw):
		bg2bw_cmd = ("bedGraphToBigWig %s %s %s"
			% (samp_bg, chrom_sizes, samp_bw))
		cmd(bg2bw_cmd, verbose)
		bigwig_err = ("Error: %s bigwig file was not generated.\n" % samp)
		if isempty(samp_bw):
			sys.exit(bigwig_err)

	# 18. call peaks
	add_dir(macs_dir)
	macs2_out=("%s/%s_peaks.broadPeak" % (macs_dir, samp))
	read_size=readsize_dic[samp]
	macs2_cmd = (("macs2 callpeak " +
		"-t %s -f BAM --keep-dup auto --nomodel " +
		"--extsize %s --broad --nolambda " +
		"-g 101274395 -n %s --outdir %s")
		% (markdups_sorted_f, read_size, samp, macs_dir))
	cmd(macs2_cmd, verbose, macs2_out)
	macs2_err = ("Error: %s did not export peaks with MACS2\n" % samp)
	if isempty(macs2_out):
		sys.exit(macs2_err)

	# 19. filter for peaks with summit <=10^-10
	add_dir(qval_dir)
	qval_f=("%s/%s_peaks.broadPeak" % (qval_dir, samp))
	broadPeak_cols=["chrom", "start", "stop", "name", "score", "strand",
		"signalVal", "pval", "qval"]
	broadPeak_df = pd.read_csv(macs2_out, sep="\t",names=broadPeak_cols)
	qval_df = broadPeak_df.loc[(
		(broadPeak_df["qval"] >= qval) & 
		(broadPeak_df["chrom"] != "ChrC") &
		(broadPeak_df["chrom"] != "ChrM")),:].copy()
	qval_df.to_csv(qval_f, header=False, index=False, sep="\t")
	broadPeak_l.append(qval_df)

# 20. concatonate all the peaks called by all the inputs
concat_df = pd.concat(broadPeak_l, ignore_index=True)
# * use the name column to group the regions called by the sample
#   from which it was called from
concat_df["name"] = (
	concat_df["name"].apply(lambda x: "_".join(x.split("_")[:-1])))
# * sort by location
concat_df.sort_values(by=["chrom","start"],
	inplace=True, ignore_index=True)
concat_df.to_csv(concat_file, header=False, index=False, sep="\t")



# 21. merge all overlapping regions and
# regions within ${merge_distance} bp apart
bedtools_merge_cmd = (("bedtools merge -c 4,5,6,7,8,9,4 " +
    "-o distinct,max,distinct,max,max,max,count_distinct " +
    "-i %s -d %i > %s")
    % (concat_file, merge_distance, merge_file))
cmd(bedtools_merge_cmd, verbose)

# 22. filter out regions that are called in less
# than ${distinct_ninputs} distinct samples
mergedRegions_df = pd.read_csv(merge_file, 
	sep="\t", names=broadPeak_cols+["ninputs"])
greenscreen_df = mergedRegions_df.loc[(
	mergedRegions_df["ninputs"] >= distinct_ninputs), :].copy()
greenscreen_df.drop(['ninputs'], axis=1, inplace=True)
greenscreen_df.to_csv(final_greenscreen, header=False, index=False, sep="\t")
if verbose:
	sys.stdout.write("> Greenscreen regions path: %s\n" % final_greenscreen)
