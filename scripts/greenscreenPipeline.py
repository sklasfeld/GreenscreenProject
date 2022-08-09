#!/usr/bin/env python

# import libraries
import os # operating system library
import sys # system-specifc library
import subprocess # subprocess management library
import tempfile # library to handle temp files
import pandas as pd # dataframe modules
import numpy as np # math modules
import argparse # argument parser library

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
def cmd(foo, quiet=False, output_f=""):
	run_cmd=True
	if len(output_f) > 0:
		if not isempty(output_f):
			run_cmd=False
	if run_cmd:
		if not quiet:
			sys.stdout.write("> CMD: %s\n" % foo)
			sys.stdout.flush()
		os.system(foo)

# Hard-Coded Variables

parser = argparse.ArgumentParser(description =
	("Generate Greenscreen Regions"))
parser.add_argument('input_meta', type = str, 
	help = ('path to file containing a comma-delimited table ' +
	'with columns for `sample_name`, `path`, and `read_length` ' +
	'for each input that will be used to make the greenscreen. ' +
	'See demo input at `demo/meta/input_meta.csv`.'))
parser.add_argument('genome_fasta', type = str, 
	help = ('path to genome sequence file in FASTA format'))
parser.add_argument('genome_gff', type = str, 
	help = ('path to genome annotation file in GFF format'))
parser.add_argument('chrom_sizes', type= str,
	help= ('path to tab-delimited file where the first column ' +
		'is the chromosome names and the second column is ' +
		'the length of the chromosome in basepairs. If this ' +
		'file does not exist yet then it will be generated.'))
parser.add_argument('out_dir', type = str, 
	help = ('path to file output directory'))
trim_arg = parser.add_argument_group('trimming arguments', 
	'arguments for trimmomatic')
trim_arg.add_argument('-ta', '--adapter_f', 
	help = ('File name containing adapters for trimming FASTQ ' +
		'reads. By default, the adapters from TRIMMOMATIC are ' +
		'used'))
trim_arg.add_argument('-tm', '--seed_mismatches', type=int, default=2,
	help = ('Trimmomatic Illuminaclip parameter which specifies ' + 
		'the maximum mismatch count which will still allow a full' +
		'match to be performed (default: 2)'))
trim_arg.add_argument('-tp', '--palindrom_clip_threshold',
	type=int, default=30, help = ('Trimmomatic Illuminaclip ' +
		'parameter which specifies how accurate the match between ' + 
		'the two `adapter ligated`  reads must be for PE palindrome ' +
		'read alignment. (default:30)'))
trim_arg.add_argument('-ts', '--simple_clip_threshold', 
	type=int, default=10, help = ('Trimmomatic Illuminaclip ' +
		'parameter which specifies how accurate the match ' +
		'between any adapter etc. sequence must be against ' + 
		'a read. (default:10)'))
trim_arg.add_argument('-t', '--trimmomatic_params', type=str,
	default="", help= ('Any other parameters to run in the ' +
		'TRIMMOMATIC command. See the TRIMMOMATIC manual ' +
		'for information on all the parameters. For the ' +
		'demo we set this to ' + 
		'`LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`'))
mapping_arg = parser.add_argument_group('mapping arguments', 
	'parameters to handle mapped reads')
mapping_arg.add_argument('-c', '--chromosomes', nargs="+",
	help = ('Filter for only reads that map to these specific ' +
		'chromosomes. For the demo set this to `Chr2`. This is ' +
		'useful if you want to ignore reads that map to the ' +
		'mitochondria, chloroplast or random scaffolds.'))
mapping_arg.add_argument('-p', '--genome_prefix', required=True,
	help = ('Base-name for genome to be used for bowtie2 indexing. ' +
		'Note that this should NOT be a path.'))
mapping_arg.add_argument('-b' ,'--bowtie2_params', type=str,
	help='Any other parameters to run in the ' +
		'Bowtie2 command. See the Bowtie2 manual ' +
		'for information on all the parameters. ' +
		'This remains unset in the demo.', default="")
mapping_arg.add_argument('-F', '--samtools_flag_filter', type=int,
	help = ('After mapping input reads to the genome, ' +
		'remove reads with these samtools flag present. ' +
		'For help to decode the meaning of SAMtools flags see: ' +
		'https://broadinstitute.github.io/picard/explain-flags.html.' +
		'(DEFAULT: 772)'), default=772)
mapping_arg.add_argument('-mq', '--samtools_mapq', type=int,
	help = ('Minimum MAPQ mapping quality. MAPQ is equal to ' +
		'−10 log10 Probability {mapping position is wrong}, ' +
		'rounded to the nearest integer. (DEFAULT: 30)'), 
	default=30)
vis_arg = parser.add_argument_group('visualization arguments', 
	'parameters to process reads into bigwig format for IGV')
vis_arg.add_argument('-s', '--scale', type=float,
	help = ('Scale normalized read value by this number for ' +
		'bigwig format. (DEFAULT: 10000000)'), default=10000000)
gs_arg = parser.add_argument_group('greenscreen arguments', 
	'settings to call greenscreen regions')
gs_arg.add_argument('-gq', '--qval', type=float,
	help = ('The maximum -log10 mean q-value across all positions ' +
		'in the broad peak regions called by MACS2. (DEFAULT: 10)'), 
	default=10)
gs_arg.add_argument('-gm', '--merge_distance', type=int,
	help = ('The maximum distance between MACS2 broad peaks for ' +
		'peaks to be merged into single greenscreen regions. ' +
		'(DEFAULT: 5000)'), default=5000)
gs_arg.add_argument('-gn', '--distinct_ninputs', type=int,
	help = ('After peaks are called in each input and regions ' +
		'withing the maximum `distinct_ninputs` are merged, regions ' +
		'which are not found significant in this set minimum of inputs ' +
		'are removed from the final greenscreen region list. '+
		'(DEFAULT: half the number of inputs used)'))
install_arg = parser.add_argument_group('installation arguments', 
	'Paths to software needed to run the script')
install_arg.add_argument('--trimmomatic_install_dir', 
	default="/usr/src/Trimmomatic-0.39", 
	help=("Path to directory which contains TRIMMATIC installation. " +
		"In the docker environment it is installed at " + 
		"`/usr/src/Trimmomatic-0.39` which is what this parameter " +
		"is set to by default."))
install_arg.add_argument('--picard_path', 
	default="/usr/src/picard/build/libs", 
	help=("Path to directory which contains PICARD installation. " +
		"In the docker environment it is installed at " + 
		"`/usr/src/picard/build/libs` which is what this parameter " +
		"is set to by default."))
out_arg = parser.add_argument_group('output arguments', 
	'set these parameters to change the output paths')
out_arg.add_argument('-ofqr','--fastqc_raw_dir', type=str,
	help=("Set this to a directory to export FASTQC files " +
	"that was run on the raw FASTA reads. By default, " +
	"these files will be output in a subdirectory called " + 
	"`FASTQC/raw` in the set `out_dir`."))
out_arg.add_argument('-oft','--fastq_trim_dir', type=str,
	help=("Set this to a directory to export trimmmed " +
		"FASTQ files after Trimmomatic. By default, " +
		"these files will be output in a subdirectory called " + 
		"`fastq/trimmed` in the set `out_dir`."))
out_arg.add_argument('-ofqt','--fastqc_trim_dir', type=str,
	help=("Set this to a directory to export FASTQC files " +
	"that was run on the trimmed FASTA reads. By default, " +
	"these files will be output in a subdirectory called " + 
	"`FASTQC/trimmed` in the set `out_dir`."))
out_arg.add_argument('-obi','--bowtie2_index_dir', type=str,
	help=("Set this to output bowtie2 genome index files in a " +
	"folder seperate from that set by `out_dir`. By default, " +
	"these files will be output in a subdirectory called " + 
	"`bowtie2_genome_dir` in the set `out_dir`."))
out_arg.add_argument('-ocs','--chipqc_samplesheet_dir', type=str,
	help=("Set this to output the chipqc samplesheets in a " +
	"folder seperate from that set by `out_dir`. By default, " +
	"these files will be output in a subdirectory called " + 
	"`chipqc_samplesheets` in the set `out_dir`."))
out_arg.add_argument('-ocr','--chipqc_out_dir', type=str,
	help=("Set this to a directory to export CHIPQC output. "+
		"By default, these files will be output in a subdirectory called " + 
	"`ChIPQCreport/inputs_noMask` in the set `out_dir`."))
out_arg.add_argument('-omp','--mapped_dir', type=str,
	help=("Set this to a directory to export all mapping data including " +
		"bowtie2 and samtools output as well as generated bigwig files. "+
		"By default, these files will be output in a subdirectory called " + 
	"`mapped/input` in the set `out_dir`."))
out_arg.add_argument('-omc','--macs_dir', type=str,
	help=("Set this to a directory to export all peak and regions files " + 
		"when finding and filtering the greenscreen regions . "+
		"By default, these files will be output in a subdirectory called " + 
		"`macs2_out/inputControls` in the set `out_dir`."))
out_arg.add_argument('-ogs','--final_greenscreen', type=str,
	help=("Set this to the file-path to output the final greenscreen " + 
		"regions. By default, this file will be output to a subdirectory " + 
		"in `macs_dir` and named based on the set greenscreen parameters"))
parser.add_argument('-q', '--quiet', default=False, action='store_true',
	help=('supress standard output'))
args = parser.parse_args()


# * check user-defined parameters
if not os.path.isfile(args.input_meta):
	sys.exit(("ERROR: Cannot find `input_meta` file at path %s. " +
		"Please set the correct file-path.") % args.input_meta)
if not os.path.isfile(args.genome_fasta):
	sys.exit(("ERROR: Cannot find `genome_fasta` file at path %s. " +
		"Please set the correct file-path.") % args.genome_fasta)
if not os.path.isfile(args.genome_gff):
	sys.exit(("ERROR: Cannot find `genome_gff` file at path %s. " +
		"Please set the correct file-path.") % args.genome_gff)
if args.adapter_f:
	if not os.path.isfile(args.adapter_f):
		sys.exit(("ERROR: Cannot find `adapter_f` file at path %s. " +
			"Please set the correct file-path.") % args.adapter_f)
if not os.path.isdir(args.out_dir):
	os.mkdir(args.out_dir)
if not args.bowtie2_index_dir:
	args.bowtie2_index_dir = ("%s/bowtie2_genome_dir" % args.out_dir)
if not args.chipqc_samplesheet_dir:
	args.chipqc_samplesheet_dir = ("%s/chipqc_samplesheet_dir" % args.out_dir)
if not args.fastqc_raw_dir:
	args.fastqc_raw_dir = ("%s/FASTQC/raw" % args.out_dir)
if not args.fastq_trim_dir:
	args.fastq_trim_dir = ("%s/fastq/trimmed" % args.out_dir)
if not args.fastqc_trim_dir:
	args.fastqc_trim_dir = ("%s/FASTQC/trimmed" % args.out_dir)
if not args.chipqc_out_dir:
	args.chipqc_out_dir = ("%s/ChIPQCreport/inputs_noMask" % args.out_dir)
if not args.mapped_dir:
	args.mapped_dir = ("%s/mapped/input" % args.out_dir)
if not args.macs_dir:
	args.macs_dir = ("%s/macs2_out/inputControls" % args.out_dir)
qval_dir = ("%s/qval%i" % (args.macs_dir, args.qval))
concat_file = ("%s/concat_input_peaks.broadPeak" % qval_dir)
merge_file = ("%s/merge%ibp_demoInputs.txt" % (qval_dir,args.merge_distance))


input_meta_df = pd.read_csv(args.input_meta)
input_meta_cols=['sample_name','path', 'read_length']
for incol in input_meta_cols:
	if incol not in list(input_meta_df.columns):
		sys.exit(("ERROR: Cannot find %s column in `input_meta` " +
			"file at path %s. This file reqires the following " +
			"columns: %s") % (incol,args.input_meta, 
			",".join(input_meta_cols)))
if len(input_meta_df["sample_name"]) != len(input_meta_df["sample_name"].unique()):
	sys.exit(("ERROR: `input_meta` file at path %s must have " +
		"unique values for each sample under `sample_name`") % 
		(args.input_meta))
input_meta_df = input_meta_df.set_index('sample_name')		

input_list=list(input_meta_df.index)

if not args.distinct_ninputs:
	args.distinct_ninputs = (
		np.floor(len(input_list)/2))
if not args.final_greenscreen:
	args.final_greenscreen = ("%s/gs_merge%ibp_call%i_demoInputs.broadPeak" 
		% (qval_dir, args.merge_distance,args.distinct_ninputs))


# 1. count genome chromosome sizes
if isempty(args.chrom_sizes):
	chrom_sizes_cmd = (
		("bash scripts/measureContigLengthFromFasta.sh %s %s")
		% (args.genome_fasta,args.chrom_sizes))
	cmd(chrom_sizes_cmd, args.quiet, args.chrom_sizes)
# * check output
if isempty(args.chrom_sizes):
	sys.exit("Error: Custom script " +
		"`scripts/measureContigLengthFromFasta.sh` " +
		"did not generate table as expected\n")

# 2. build genome bowtie2 index
add_dir(args.bowtie2_index_dir)
bowtie2_index_prefix  =  ("%s/%s" % (args.bowtie2_index_dir, args.genome_prefix))
if not os.path.isdir(args.bowtie2_index_dir):
	bowtie2_index_cmd = (
		"bowtie2-build %s %s" %
		(args.genome_fasta, bowtie2_index_prefix))
	cmd(bowtie2_index_cmd, args.quiet)
elif len(os.listdir(args.bowtie2_index_dir)) == 0:
	bowtie2_index_cmd = (
		"bowtie2-build %s %s" %
		(args.genome_fasta, bowtie2_index_prefix))
	cmd(bowtie2_index_cmd, args.quiet)
bowtie2_index_err=("Error: Bowtie2 did not build genome index\n")
# * check output
if not os.path.isdir(args.bowtie2_index_dir):
	sys.exit(bowtie2_index_err)
if len(os.listdir(args.bowtie2_index_dir)) == 0:
	sys.exit(bowtie2_index_err)


round=1
ninputs=len(input_list)
broadPeak_l=[]
for samp, samp_meta_df in input_meta_df.iterrows():
	raw_fastq = samp_meta_df["path"]
	
	# 3. run fastqc on raw fastq
	add_dir(args.fastqc_raw_dir)
	fastqc_raw_output = ("%s/%s_fastqc.html" %
		(args.fastqc_raw_dir, samp))
	fastqc_raw_cmd = (("fastqc -o %s %s") % 
		(args.fastqc_raw_dir, raw_fastq))
	cmd(fastqc_raw_cmd, args.quiet, fastqc_raw_output)
	# * check output
	
	if isempty(fastqc_raw_output):
		sys.exit("Error: FASTQC output for %s reads were not generated\n"
			% samp)
	
	# 4. run trimming
	copied_adapter=False
	new_adapter_dir=False
	temp_adapter_dir="adapters"
	if not args.adapter_f:
		copied_adapter=True
		add_dir(temp_adapter_dir)
		if len(os.listdir(temp_adapter_dir)) == 0:
			new_adapter_dir=True
		args.adapter_f = ("%s/TruSeq3-SE.fa" % temp_adapter_dir)
		cp_adapters_cmd = ("cp %s/adapters/TruSeq3-SE.fa %s"
			% (args.trimmomatic_install_dir, args.adapter_f))
		cmd(cp_adapters_cmd, args.quiet, args.adapter_f)
		
	
	add_dir(args.fastq_trim_dir)
	trim_fastq = ("%s/%s.trimmed.fq.gz" % (args.fastq_trim_dir,samp))
	trim_cmd = (("java -jar %s/trimmomatic-0.39.jar "+
		"SE -phred33 %s %s ILLUMINACLIP:%s:%i:%i:%i %s TOPHRED33") % 
		(args.trimmomatic_install_dir, raw_fastq, trim_fastq, 
			args.adapter_f, args.seed_mismatches, 
			args.palindrom_clip_threshold, args.simple_clip_threshold,
			args.trimmomatic_params))
	cmd(trim_cmd, args.quiet, trim_fastq)
	# * check output
	if isempty(trim_fastq):
		sys.exit("Error: Trimmomatic output for %s was not generated\n"
			% samp)

	# 5. run fastqc on trimmed fastq
	add_dir(args.fastqc_trim_dir)
	fastqc_trim_output = ("%s/%s.trimmed_fastqc.html" %
		(args.fastqc_trim_dir, samp))
	fastqc_trim_cmd = ("fastqc -o %s %s"
		%(args.fastqc_trim_dir, trim_fastq))
	cmd(fastqc_trim_cmd, args.quiet, fastqc_trim_output)
	# * check output
	if isempty(fastqc_trim_output):
		sys.exit("Error: FASTQC output for %s trimmed reads were not generated\n"
			% samp)
	
	# 6. map reads to genome with bowtie2
	add_dir(args.mapped_dir)
	bowtie2_mapped_f = ("%s/%s.sam" %
		(args.mapped_dir, samp))
	bowtie2_map_cmd = (("bowtie2 %s --phred33 -q -x %s -S %s %s") % 
		(args.bowtie2_params, bowtie2_index_prefix, bowtie2_mapped_f, trim_fastq))
	cmd(bowtie2_map_cmd, args.quiet, bowtie2_mapped_f)
	# * check output
	bowtie2_mapped_file_err = (("Error: bowtie2 mapping output for %s was " +
		"not generated\n") % samp)
	if isempty(bowtie2_mapped_f):
		sys.exit(bowtie2_mapped_file_err)


	# 7. sort the reads
	bowtie2_mapped_sorted_f = ("%s/%s.bam" % (args.mapped_dir, samp))
	bowtie2_mapped_sorted_cmd = (
		"samtools sort -o %s %s" % 
		(bowtie2_mapped_sorted_f, bowtie2_mapped_f))
	cmd(bowtie2_mapped_sorted_cmd, args.quiet, bowtie2_mapped_sorted_f)
	# * check output
	bowtie2_mapped_sorted_file_err = (("Error: bowtie2 mapping output for " +
		"%s was not sorted\n") % samp)
	if isempty(bowtie2_mapped_sorted_f):
		sys.exit(bowtie2_mapped_sorted_file_err)

	# 8. index the bam file
	bowtie2_mapped_sorted_idx = ("%s.bai" % bowtie2_mapped_sorted_f)
	index_sorted_bowtie2_mapped_cmd = (
		"samtools index %s" % bowtie2_mapped_sorted_f)
	cmd(index_sorted_bowtie2_mapped_cmd, args.quiet, 
		bowtie2_mapped_sorted_idx)
	# * check output
	bowtie2_mapped_sorted_idx_err = (("Error: sorted bowtie2 mapping " +
		"output for %s was not indexed\n") % samp)
	if isempty(bowtie2_mapped_sorted_idx):
		sys.exit(bowtie2_mapped_sorted_idx_err)

	# 9. remove reads without MAPQ>=30 and without sam flags 772
	bowtie2_filtered_f = ("%s/%s.filter.bam" % (args.mapped_dir, samp))
	bowtie2_samtools_filter_cmd = (("samtools view -F %i " + 
		"-q %i -b %s -o %s") % (args.samtools_flag_filter, 
		args.samtools_mapq, bowtie2_mapped_sorted_f, bowtie2_filtered_f))
	if args.chromosomes:
		bowtie2_samtools_filter_cmd = ("%s %s" % 
			(bowtie2_samtools_filter_cmd, "_".join(args.chromosomes)))
	cmd(bowtie2_samtools_filter_cmd, args.quiet, bowtie2_filtered_f)
	# * check output
	bowtie2_filter_file_err = (("Error: bowtie2 mapped output for " +
		"%s could not be filtered (MAPQ >= %i and only include " +
		"reads with none of the FLAGS in %i present)\n") % 
		(samp, args.samtools_mapq, args.samtools_flag_filter))
	if isempty(bowtie2_filtered_f):
		sys.exit(bowtie2_filter_file_err)

	# 10. sort filtered reads
	bowtie2_sorted_filtered_f = ("%s/%s.filter.sorted.bam" % (args.mapped_dir, samp))
	bowtie2_sort_filtered_cmd = ("samtools sort %s -o %s" 
		% (bowtie2_filtered_f, bowtie2_sorted_filtered_f))
	cmd(bowtie2_sort_filtered_cmd, args.quiet,
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
	cmd(index_sorted_filtered_bowtie2_cmd, args.quiet, 
		bowtie2_sorted_filtered_idx)
	# * check output
	bowtie2_sorted_filtered_idx_err = (("Error: sorted filtered mapped reads " +
		"from bowtie2 for %s was not indexed\n") % samp)
	if isempty(bowtie2_sorted_filtered_idx):
		sys.exit(bowtie2_sorted_filtered_idx_err)

	# 12. mark duplicates with picard
	picard_out = ("%s/%s.dupmark.bam" % (args.mapped_dir, samp))
	markdups_cmd = (("java -jar %s/picard.jar MarkDuplicates " +
		"-I %s -O %s " +
		"-M %s/%s.dup.qc " +
		"-VALIDATION_STRINGENCY LENIENT " +
		"-REMOVE_DUPLICATES false " +
		"-ASSUME_SORTED true") %
		(args.picard_path, bowtie2_sorted_filtered_f, picard_out,
		args.mapped_dir, samp))
	cmd(markdups_cmd, args.quiet, picard_out)
	# * check output
	markdups_err = (("Error: duplicate were not marked for %s\n") % samp)
	if isempty(picard_out):
		sys.exit(markdups_err)

	# 13. sort reads after marking the duplicates
	markdups_sorted_f = ("%s/%s.dupmark.sorted.bam" % (args.mapped_dir, samp))
	markdups_sort_cmd = ("samtools sort %s -o %s" 
		% (picard_out, markdups_sorted_f))
	cmd(markdups_sort_cmd, args.quiet, markdups_sorted_f)
	# * check output
	markdups_sorted_file_err = (("Error: picard markdups output for " +
		"%s was not sorted\n") % samp)
	if isempty(markdups_sorted_f):
		sys.exit(markdups_sorted_file_err)


	# 14. Generate index for the sorted marked reads
	markdups_sorted_idx = ("%s.bai" % markdups_sorted_f)
	index_markdups_cmd = (
		"samtools index %s" % markdups_sorted_f)
	cmd(index_markdups_cmd, args.quiet, markdups_sorted_idx)
	# * check output
	markdups_idx_err = ("Error: picard output for %s was not indexed\n" % samp)
	if isempty(markdups_sorted_idx):
		sys.exit(markdups_idx_err)


	# 15. Generate sample sheet for ChIPQC code
	add_dir(args.chipqc_samplesheet_dir)
	samplesheet=("%s/chipqc_sampleSheet_%s.csv" % (args.chipqc_samplesheet_dir, samp))
	samplesheet_handle = open(samplesheet, "w")
	samplesheet_handle.write("SampleID,bamReads\n")
	samplesheet_handle.write("%s,%s/%s.dupmark.sorted.bam\n" % 
		(samp, args.mapped_dir,samp))
	samplesheet_handle.close()

	# 15. Run ChIPQC
	add_dir(args.chipqc_out_dir)
	chipqc_out = (("%s/%s/ChIPQC.html") % (args.chipqc_out_dir, samp))
	chipqc_r_cmd = "Rscript scripts/ChIPQC.R --indivReports -g Araport11 "
	if args.chromosomes:
		chipqc_r_cmd = ("%s -c %s" % (chipqc_r_cmd, "_".join(args.chromosomes)))
	chipqc_r_cmd = (("%s -a %s -s %s %s %s") % 
		(chipqc_r_cmd, args.genome_gff, 
			args.chrom_sizes, samplesheet, args.chipqc_out_dir))
	cmd(chipqc_r_cmd, args.quiet, chipqc_out)

	# 16. normalize signal and output BEDGRAPH
	totreads= bamsize(markdups_sorted_f)
	scaling=args.scale / totreads
	
	samp_bg = ("%s/%s.bg" % (args.mapped_dir, samp))
	if isempty(samp_bg):
		temp_bg = tempfile.NamedTemporaryFile()
		scale_bam = ("genomeCoverageBed -ibam %s -bg -scale %f > %s" 
			% (markdups_sorted_f, scaling, temp_bg.name))
		cmd(scale_bam, args.quiet)
	
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
	samp_bw = ("%s/%s.bw" % (args.mapped_dir, samp))
	if not os.path.isfile(samp_bw):
		bg2bw_cmd = ("bedGraphToBigWig %s %s %s"
			% (samp_bg, args.chrom_sizes, samp_bw))
		cmd(bg2bw_cmd, args.quiet)
		bigwig_err = ("Error: %s bigwig file was not generated.\n" % samp)
		if isempty(samp_bw):
			sys.exit(bigwig_err)

	# 18. call peaks
	add_dir(args.macs_dir)
	macs2_out=("%s/%s_peaks.broadPeak" % (args.macs_dir, samp))
	read_size = input_meta_df.loc[samp,"read_length"]
	macs2_cmd = (("macs2 callpeak " +
		"-t %s -f BAM --keep-dup auto --nomodel " +
		"--extsize %s --broad --nolambda " +
		"-g 101274395 -n %s --outdir %s")
		% (markdups_sorted_f, read_size, samp, args.macs_dir))
	cmd(macs2_cmd, args.quiet, macs2_out)
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
		(broadPeak_df["qval"] >= args.qval) & 
		(broadPeak_df["chrom"] != "ChrC") &
		(broadPeak_df["chrom"] != "ChrM")),:].copy()
	qval_df.to_csv(qval_f, header=False, index=False, sep="\t")
	broadPeak_l.append(qval_df)

# remove temporary adapter file/directory
if copied_adapter:
	# remove adapter file
	if not args.quiet:
		print("remove %s" % args.adapter_f)
	os.remove(args.adapter_f)
	if new_adapter_dir:
		if not args.quiet:
			print("remove %s" % temp_adapter_dir)
		os.rmdir(temp_adapter_dir)
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
    % (concat_file, args.merge_distance, merge_file))
cmd(bedtools_merge_cmd, args.quiet)

# 22. filter out regions that are called in less
# than ${distinct_ninputs} distinct samples
mergedRegions_df = pd.read_csv(merge_file, 
	sep="\t", names=broadPeak_cols+["ninputs"])
greenscreen_df = mergedRegions_df.loc[(
	mergedRegions_df["ninputs"] >= args.distinct_ninputs), :].copy()
greenscreen_df.drop(['ninputs'], axis=1, inplace=True)
greenscreen_df.to_csv(args.final_greenscreen, header=False, index=False, sep="\t")
if not args.quiet:
	sys.stdout.write("> Greenscreen regions path: %s\n" % args.final_greenscreen)
