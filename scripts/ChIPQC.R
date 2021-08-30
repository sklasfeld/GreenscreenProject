
#!/usr/bin/env Rscript

# 2021, Samantha Klasfeld, the Wagner Lab
# the Perelman School of Medicine, the University of Pennsylvania
# Samantha Klasfeld, 7-28-2021

library("argparse")
options(error=traceback)
argv = commandArgs(trailingOnly=TRUE)

runchipqc <- function(exp_obj, out_path, facet=NULL, colour=NULL){
  if ((!is.null(facet)) &&
      (!is.null(colour))){
    ChIPQCreport(exp_obj,facetBy=facet,
                 colourBy = colour,
                 reportFolder=out_path,
                 reportName=paste0("ChIPQC"))
  } else if (!is.null(facet)){
    ChIPQCreport(exp_obj,facetBy=facet,
                 reportFolder=out_path,
                 reportName=paste0("ChIPQC"))
  } else if (!is.null(colour)){
    ChIPQCreport(exp_obj,
                 colourBy = colour,
                 reportFolder=out_path,
                 reportName=paste0("ChIPQC"))
  } else{
    ChIPQCreport(exp_obj,
                 reportFolder=out_path,
                 reportName=paste0("ChIPQC"))
  }
}

parser <- ArgumentParser(description=paste('Create ChIPQC reports',
  'for each sample in a sample sheet'))
parser$add_argument('exp_design',  help=paste('A file containing',
  'a csv table describing the ChIP experiment samples. Required',
  'column names include: SampleID, bamReads, and Peaks. Other',
  'column should futher describe the experient (ie. replicate',
  'number, treatment, genotype).'))
parser$add_argument('out_dir',  help='directory for output')
parser$add_argument('-g','--genome', help = paste('Genome',
  'name. The following annotation specifiers are supported:',
  '"hg19" (Human, version 19), "hg18" (Human, version 18),',
  '"mm10" (Mouse, version 10), "mm9" (Mouse, version 9),',
  '"rn4" (Rat, version 4), "ce6" (C. Elgans, version 6),',
  '"dm3" (D. Melanogaster, version 3). For other genomes',
  '(eg.Araport11) one must set a custom annotation',
  'using the `gff_ann` parameter. (default:hg19)'), default="hg19")
parser$add_argument('-c','--chrom', nargs="+", 
  help = paste('Specification of which chromosomes to use for',
  'computing QC statistics. If missing, the first chromosome',
  'which has a peak is checked. If NULL, all chromosomes will be',
  'checked (which may be time-consuming).'))
parser$add_argument('-a','--gff_ann',  help=paste('set annotation',
  'file in GFF3 format if using custom genome annotation'))
parser$add_argument('-s','--chr_sizes', 
  help=paste('file containing a tab-delimited table. The',
  'first column contains the chromosome/scaffold names.',
  'The second column contains the lengths of the respective',
  'chromosome/scaffolds. This is required if using custom genome',
  'annotation.'))
parser$add_argument('-f', '--facet', 
  help = "column to group the samples by")
parser$add_argument('-k','--color', 
  help = "column to color the samples by")
parser$add_argument("-i", "--indivReport#!/bin/bash

# Goal: scale-down a bam file to a specific number of reads
# through random subsampling of the total reads

if [ $# -lt 3 ]; then
   echo 'Command Format: subsample_bam_V1.sh <inbam> <n> <outbam>'
   exit 1
fi


# seed (this is set for reproducibility purposes)
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

# import the bam file and the number the scale down to
inbam=${1%.bam}
n=$2
# import the name of the output bam file
outbam=${3%.bam}

# check that the inbam file exists
if [ ! -f "${inbam}.bam" ]; then
	echo "ERROR: ${inbam}.bam does not exist"
	exit 1
fi
# check that the outbam.sam file does not exist
if [ -f "${outbam}.sam" ]; then
	echo "WARNING: Overwritting ${outbam}.sam"
fi
if [ "${inbam}.bam" == "${outbam}.bam" ]; then
	echo "ERROR: INPUT AND OUTPUT FILE CANNOT BE THE SAME"
	exit 1
fi

# get header for output sam file
samtools view -H ${inbam}.bam > ${outbam}.sam
# 1. convert bam file to human readable sam format
# 2. shuffle the reads
# 3. return the first  $n reads
samtools view ${inbam}.bam | shuf --random-source=<(get_seeded_random 42) - | head -n ${n} >> ${outbam}.sam

# convert the sam file to bam format
samtools sort -O BAM ${outbam}.sam > ${outbam}.bam
samtools index ${outbam}.bam
s", action="store_true",
  default=F, help="Print a seperate report for each sample")
args <- parser$parse_args()

require(ChIPQC)
require(GenomicFeatures)
require(GenomicRanges)


# check if need custom gff_ann file
supported_genomes=c("hg19","hg18","mm10","mm9","rn4","ce6","dm3")
if (!(args$genome %in% supported_genomes) && 
    (is.null(args$gff_ann) || is.null(args$chr_sizes))){
  stop("`gff_ann` and `chr_sizes` parameters must be set for custom genomes")
} else{
    if (!(is.null(args$chr_sizes))){
        chrSize_table = read.table(file.path(args$chr_sizes), sep="\t", header=F)
        colnames(chrSize_table)=c("seqnames","seqlengths")
        chrSize_seqInfo = Seqinfo(seqnames=chrSize_table$seqnames,
            seqlengths=chrSize_table$seqlengths,
            isCircular=rep(c(FALSE), length(chrSize_table$seqnames)),
            genome=args$genome)
    }
    samples = read.csv(file.path(args$exp_design), sep=",")


    # check `exp_design` table has correct column names
    if (!("SampleID" %in% colnames(samples) && 
          "bamReads" %in% colnames(samples))){
        stop(paste('`exp_design` table must contain the columns:',
            '`SampleID`, `Peaks`, and `bamReads`'))
    } else {
        missing_files_str=c()
        no_peaks=c()
        if (!("Peaks" %in% colnames(samples))){
            samples[,"Peaks"] <- NA
        } else{
            for(row in 1:nrow(samples)){
                if (!file.exists(samples$Peaks[row]) &&
                    length(samples$Peaks[row]) > 0){
                    missing_files_str = \
                      append(missing_files_str,samples$Peaks[row])
                }else{
                    if(file.info(samples$Peaks[row])$size == 0){
                        no_peaks = append(no_peaks,samples$Peaks[row])
                    }
                }
            }
        }
        for(row in 1:nrow(samples)){
            if (!file.exists(samples$bamReads[row])){
                missing_files_str = \
                  append(missing_files_str,samples$bamReads[row])
            }
            index_file = paste0(samples$bamReads[row],".bai")
            if (!file.exists(index_file)){
                missing_files_str = \
                  append(missing_files_str,index_file)
            }
        }
        # check files input into script
        if (length(missing_files_str)>0){
            stop(paste0("Error in exp_design file (",
                args$exp_design,"). Cannot find the ",
                "following files:"))
            for (f in missing_files_str){
                print(f)
            }
        } else{
            if (!(is.null(args$gff_ann))){
                # GIVEN A GENOME WITH CUSTOM ANNOTATION
                txdb <- makeTxDbFromGFF(args$gff_ann, 
                                        format="gff3",
                                        chrominfo = chrSize_seqInfo) 
                All5utrs <- (
                  reduce(unique(unlist(fiveUTRsByTranscript(txdb)))))
                All3utrs <- (
                  reduce(unique(unlist(threeUTRsByTranscript(txdb)))))
                Allcds <- (
                  reduce(unique(unlist(cdsBy(txdb,"tx")))))
                Allintrons <- (
                  reduce(unique(unlist(intronsByTranscript(txdb)))))
                Alltranscripts <- transcripts(txdb)
                
                posAllTranscripts <- (
                  Alltranscripts[strand(Alltranscripts) == "+"])
                posAllTranscripts <- (
                  posAllTranscripts[!(start(posAllTranscripts)-20000 
                  < 0)])
                negAllTranscripts <- (
                  Alltranscripts[strand(Alltranscripts) == "-"])
                chrLimits <- (
                 seqlengths(negAllTranscripts)[as.character(
                 seqnames(negAllTranscripts))])    
                negAllTranscripts <- (
                  negAllTranscripts[!(end(negAllTranscripts)+20000 > 
                  chrLimits)])      
                Alltranscripts <- c(posAllTranscripts,negAllTranscripts)
                Promoters500 <-  reduce(flank(Alltranscripts,500))    
                Promoters2000to500 <-  reduce(flank(Promoters500,1500))
                LongPromoter20000to2000  <- (
                  reduce(flank(Promoters2000to500,18000)))
                
                customAnnotation <- list(version=args$genome,
                  LongPromoter20000to2000=LongPromoter20000to2000,
                  Promoters2000to500=Promoters2000to500,
                  Promoters500=Promoters500, All5utrs=All5utrs,
                  Alltranscripts=Alltranscripts, Allcds=Allcds,
                  Allintrons=Allintrons,All3utrs=All3utrs)
                
                # remove all peak files that do not contain
                # peaks from the samples table
                if (length(no_peaks) > 0){
                    samples[which(
                      samples$Peaks %in% no_peaks),"Peaks"]=NA
                }
                
                if (!(is.null(args$chrom))){
                    exp_obj=ChIPQC(samples, customAnnotation,
                        chromosomes = args$chrom)
                } else {
                    exp_obj=ChIPQC(samples, customAnnotation,
                        chromosome=NULL)
                }
                out_path=file.path(paste0(args$out_dir))
                
                if(! args$indivReports){
                  print(">> write main report")
                  runchipqc(exp_obj, out_path, args$facet, args$color)
                } else{
                  print(">> write sample reports")
                  for (row in 1:nrow(samples)) {
                    print(paste(">>",samples[row]))
                    if (!(is.null(args$chrom))){
                      exp_obj_#!/bin/bash

working_dir="mapped/chip"
genome="meta/Araport11/TAIR10_chr_count.txt"

# return average signal of two reps
pool_two=("FD_W_2020" "TFL1_A_W_2020"
  "TFL1_B_W_2020" "TFL1_fd_W_2020"
  "LFY_P_2016" "FD_S_2019" 
  "FD_ft10_tsf1_S_2019" "LFY_P_2011"
  "FD_C_2020")
for samp in "${pool_two[@]}"; do
  bedtools unionbedg -i \
    ${working_dir}/${samp}_R1.bg \
    ${working_dir}/${samp}_R2.bg \
    > ${working_dir}/${samp}.unionbg
  awk 'BEGIN{OFS="\t"} \
    $1!="ChrC" && $1!="ChrM"{ \
    avg=($4+$5)/2;
    print $1,$2,$3,avg}' \
    ${working_dir}/${samp}.unionbg > \
    ${working_dir}/${samp}.bedgraph
  bedGraphToBigWig ${working_dir}/${samp}.bedgraph \
    ${genome} ${working_dir}/${samp}.bigwig
done

# return average signal of three reps
pool_three=("TFL1_S_2020" "LFY_W_2021")
for samp in "${pool_three[@]}"; do
  bedtools unionbedg -i \
    ${working_dir}/${samp}_R1.bg \
    ${working_dir}/${samp}_R2.bg \
    ${working_dir}/${samp}_R3.bg \
    > ${working_dir}/${samp}.unionbg
  awk 'BEGIN{OFS="\t"} \
    $1!="ChrC" && $1!="ChrM"{ \
    avg=($4+$5+$6)/3;
    print $1,$2,$3,avg}' \
    ${working_dir}/${samp}.unionbg > \
    ${working_dir}/${samp}.bedgraph
  bedGraphToBigWig ${working_dir}/${samp}.bedgraph \
    ${genome} ${working_dir}/${samp}.bigwig
donei = ChIPQCsample(
                        reads=file.path(samples[row,"bamReads"]),
                        peaks = file.path(samples[row,"Peaks"]),
                        annotation = customAnnotation,
                        runCrossCor=T, chromosomes = args$chrom)
                    } else{
                      exp_obj_i = ChIPQCsample(
                        reads=file.path(samples[row,"bamReads"]),
                        peaks = file.path(samples[row,"Peaks"]),
                        annotation = customAnnotation,
                        runCrossCor=T, chromosomes=NULL)
                    }
                    out_path = file.path(
                        paste0(args$out_dir,"/",samples[row,"SampleID"]))
                    runchipqc(exp_obj, out_path, 
                        args$facet, args$color)
                  }
                }
            } else{
                # GIVEN A GENOME WITH DEFAULT ANNOTATION
                if (length(no_peaks) > 0){
                    samples[which(samples$Peaks %in% 
                    no_peaks),"Peaks"]=NA
                }
                
                if (!(is.null(args$chrom))){
                    exp_obj=ChIPQC(samples, args$genome,
                        chromosomes = args$chrom)
                } else {
                    exp_obj=ChIPQC(samples, args$genome)
                }
                out_path=file.path(paste0(args$out_dir))
                
                if(! args$indivReports){
                    print(">> write main report")
                    runchipqc(exp_obj, out_path, 
                        args$facet, args$color)
                } else{
                    print(">> write sample reports")
                    for (row in 1:nrow(samples)) {
                        print(paste(">>",samples[row]))
                        if (!(is.null(args$chrom))){
                            exp_obj_i = ChIPQCsample(
                                reads=file.path(samples[row,"bamReads"]),
                                peaks = file.path(samples[row,"Peaks"]),
                                annotation = args$genome,
                                runCrossCor=T, chromosomes = args$chrom)
                        } else{
                            exp_obj_i=ChIPQCsample(
                                reads=file.path(samples[row,"bamReads"]),
                                peaks = file.path(samples[row,"Peaks"]),
                                annotation = args$genome,
                                runCrossCor=T, chromosomes=NULL)
                        }
                        out_path = file.path(
                            paste0(args$out_dir,"/",samples[row,"SampleID"]))
                        runchipqc(exp_obj, out_path, 
                            args$facet, args$color)
                    }
                }
            }
        }
    }
}
