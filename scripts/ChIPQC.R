#!/usr/bin/env Rscript

# 2021, Samantha Klasfeld, the Wagner Lab
# the Perelman School of Medicine, the University of Pennsylvania
# Samantha Klasfeld, 7-28-2021

# install argparse
if (!("argparse" %in% rownames(installed.packages()))){
    install.packages("argparse", repos='http://cran.us.r-project.org')
}


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

  if (!is.null(colour)){
    #update_geom_defaults("line", list(size = .5, colour=c('black','darkgreen','blue')))
    update_geom_defaults("line", list(size = .5, alpha=.5))
    q <- plotCC(exp_obj,
      colourBy=colour, facet=FALSE) +
      scale_color_manual(values=c('black','darkgreen','blue')) 
    ggsave(plot = q, width = 4, height = 4, 
      filename = file.path(out_path,"CCplot2.png"))
  }
}

parser <- ArgumentParser(description=paste('Create ChIPQC reports for each',
  'sample in a sample sheet'))
parser$add_argument('exp_design',  help=paste('A file containing a csv table',
  'describing the ChIP experiment samples. Required column names include:',
  'SampleID, bamReads, and Peaks. Other column should futher describe the',
  'experient (ie. replicate number, treatment, genotype).'))
parser$add_argument('out_dir',  help='directory for output')
parser$add_argument('-g','--genome', help = paste('Genome name. The following',
    'annotation specifiers are supported: "hg19" (Human, version 19),',
    '"hg18" (Human, version 18), "mm10" (Mouse, version 10),',
    '"mm9" (Mouse, version 9),  "rn4" (Rat, version 4),',
    '"ce6" (C. Elgans, version 6),  "dm3" (D. Melanogaster, version 3).',
    'For other genomes (eg.Araport11) one must set a custom annotation',
    'using the `gff_ann` parameter. (default:hg19)'), default="hg19")
parser$add_argument('-c','--chrom', nargs="+", help = paste('Specification of',
  'which chromosomes to use for computing QC statistics. If missing, the first',
  'chromosome which has a peak is checked. If NULL, all chromosomes will be',
  'checked (which may be time-consuming).'))
parser$add_argument('-a','--gff_ann',  help=paste('set annotation file in GFF3',
    'format if using custom genome annotation'))
parser$add_argument('-s','--chr_sizes', help=paste('file containing a',
    'tab-delimited table. The first column contains the chromosome/scaffold',
    'names.The second column contains the lengths of the respective',
    'chromosome/scaffolds. This is required if using custom genome',
    'annotation.'))
parser$add_argument('-k','--color', help = paste("column to color the samples by"))
parser$add_argument("-i", "--indivReports", dest="indivReports",
                     action="store_true", default=F,
                    help="Print a seperate report for each row in exp_design sheet")
parser$add_argument("-r", "--sampleReports", dest="sampleReports",
                    action="store_true", default=F,
                    help=paste('Print a seperate report for each facet value (set by',
                        'the `facet` parameter) in the exp_design sheet.',
                        'This cannot be set with `indivReports` parameter.'))
parser$add_argument('-f', '--facet', help = paste("column to group the samples by"), 
    default='SampleID')
args <- parser$parse_args()

# install CRAN packages
if( !("rsvg" %in% rownames(installed.packages())) ) {
    install.packages("rsvg",
        repos='https://cran.r-project.org')
}

# get Bioconductor installation if not available
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", dependencies=TRUE,
        repos='http://cran.us.r-project.org')
}


# get Bioconductor libraries installed if not available
listOfBiocPackages=c('ShortRead', 'edgeR', 'DESeq2', 'GOstats',
    'amap', 'systemPipeR', 'ChIPQC',
    'GenomicFeatures','GenomicRanges')

notInstalledBiocPackages <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
## check there's still something left to install
if( length(notInstalledBiocPackages) ) {
    for (biocPackage in listOfBiocPackages[ notInstalledBiocPackages ]) {
        BiocManager::install(biocPackage)
    }
}

print("> All packages should be installed")

require(ChIPQC)
require(GenomicFeatures)
require(GenomicRanges)

print("> All packages are loaded.")

print("> Preparing objects for ChIPQC")

# check that both params are not set
if ((args$indivReports == TRUE) && (args$sampleReports == TRUE)){
    stop('params `indivReports` and `sampleReports` cannot both be set')
} else{
    # check if need custom gff_ann file
    supported_genomes=c("hg19","hg18","mm10","mm9","rn4","ce6","dm3")
    if (!(args$genome %in% supported_genomes) && 
        (is.null(args$gff_ann) || is.null(args$chr_sizes))){
      stop("`gff_ann` and `chr_sizes` parameters must be set for custom genomes")
    } else{
        if (!(is.null(args$chr_sizes))){
            chrSize_table = read.table(file.path(args$chr_sizes), sep="\t", header=F)
            colnames(chrSize_table)=c("seqnames","seqlengths")
            chrSize_seqInfo = Seqinfo(
                seqnames=as.character(chrSize_table$seqnames),
                seqlengths=chrSize_table$seqlengths,
                isCircular=rep(c(FALSE), length(chrSize_table$seqnames)),
                genome=args$genome)
        }
        samples = read.csv(file.path(args$exp_design), sep=",")


        # check `exp_design` table has correct column names
        if (!("SampleID" %in% colnames(samples) && 
              "bamReads" %in% colnames(samples))){
            stop(paste('`exp_design` table must contain the columns:',
                '`SampleID` and `bamReads`'))
        } else {
            missing_files_str=c()
            no_peaks=c()
            noPeakFiles=FALSE
            if (!("Peaks" %in% colnames(samples))){
                samples[,"Peaks"] <- NA
                noPeakFiles=TRUE
            } else{
                for(row in 1:nrow(samples)){
                    peaks_file=file.path(samples$Peaks[row])
                    if (!file.exists(peaks_file) &&
                        length(samples$Peaks) > 0){
                        missing_files_str = append(missing_files_str,
                            file.path(samples$Peaks[row]))
                    }else{
                        if(file.info(samples$Peaks[row])$size == 0){
                            no_peaks = append(no_peaks,samples$Peaks[row])
                        }
                    }
                }
            }
            for(row in 1:nrow(samples)){
                if (!file.exists(file.path(samples$bamReads[row]))){
                    missing_files_str = append(missing_files_str,
                        file.path(samples$bamReads[row]))
                }
                index_file = paste0(
                    file.path(samples$bamReads[row]),".bai")
                if (!file.exists(index_file)){
                    missing_files_str = append(missing_files_str,index_file)
                }
            }
            # check files input into script
            if (length(missing_files_str)>0){
                err_str=paste0("Error in exp_design file (",
                    args$exp_design,"). Cannot find the ",
                    "following files: ")
                for (f in missing_files_str){
                    err_str=paste0(err_str,f,", ")
                }
                stop(err_str)
            } else{
                if (!(is.null(args$gff_ann))){
                    # GIVEN A GENOME WITH CUSTOM ANNOTATION
                    txdb <- makeTxDbFromGFF(args$gff_ann, 
                                            format="gff3",
                                            chrominfo = chrSize_seqInfo) 
                    All5utrs <- reduce(unique(unlist(fiveUTRsByTranscript(txdb))))
                    All3utrs <- reduce(unique(unlist(threeUTRsByTranscript(txdb))))
                    Allcds <- reduce(unique(unlist(cdsBy(txdb,"tx"))))
                    Allintrons <- reduce(unique(unlist(intronsByTranscript(txdb))))
                    Alltranscripts <- transcripts(txdb)
                    
                    posAllTranscripts <- Alltranscripts[strand(Alltranscripts) == "+"]
                    posAllTranscripts <- posAllTranscripts[!(start(posAllTranscripts)-20000 < 0)]
                    negAllTranscripts <- Alltranscripts[strand(Alltranscripts) == "-"]
                    chrLimits <- seqlengths(negAllTranscripts)[as.character(seqnames(negAllTranscripts))]      
                    negAllTranscripts <- negAllTranscripts[!(end(negAllTranscripts)+20000 > chrLimits)]      
                    Alltranscripts <- c(posAllTranscripts,negAllTranscripts)
                    Promoters500 <-  reduce(flank(Alltranscripts,500))    
                    Promoters2000to500 <-  reduce(flank(Promoters500,1500))
                    LongPromoter20000to2000  <- reduce(flank(Promoters2000to500,18000))
                    
                    customAnnotation <- list(version=args$genome,
                      LongPromoter20000to2000=LongPromoter20000to2000,
                      Promoters2000to500=Promoters2000to500,Promoters500=Promoters500,
                      All5utrs=All5utrs,Alltranscripts=Alltranscripts,Allcds=Allcds,
                      Allintrons=Allintrons,All3utrs=All3utrs)
                    
                    # remove all peak files that do not contain
                    # peaks from the samples table
                    if (length(no_peaks) > 0){
                        samples[which(samples$Peaks %in% no_peaks),"Peaks"]=NA
                    }
                    out_path=file.path(paste0(args$out_dir))
                    if((! args$indivReports) && (! args$sampleReports)){
                        print(">> write main report given custom annotation")
                        if (!(is.null(args$chrom))){
                            exp_obj=ChIPQC(samples, customAnnotation,
                            chromosomes = args$chrom)
                        } else {
                            exp_obj=ChIPQC(samples, customAnnotation,
                                chromosome=NULL)
                        }
                        runchipqc(exp_obj, out_path, args$facet, args$color)
                    } else if( args$indivReports){
                        print(">> write individual reports given custom annotation")
                        for (row in 1:nrow(samples)) {
                            print(paste(">>",samples[row,"SampleID"]))
                            if (! (noPeakFiles)){
                                if (!(is.null(args$chrom))){
                                    exp_obj_i = ChIPQCsample(
                                        reads=file.path(samples[row,"bamReads"]),
                                        peaks = file.path(samples[row,"Peaks"]),
                                        annotation = customAnnotation,
                                        runCrossCor=T, chromosomes = args$chrom)
                                } else{
                                    exp_obj_i=ChIPQCsample(
                                        reads=file.path(samples[row,"bamReads"]),
                                        peaks = file.path(samples[row,"Peaks"]),
                                        annotation = customAnnotation,
                                        runCrossCor=T, chromosomes=NULL)
                                }
                            } else{
                                if (!(is.null(args$chrom))){
                                    exp_obj_i = ChIPQCsample(
                                        reads=file.path(samples[row,"bamReads"]),
                                        annotation = customAnnotation,
                                        runCrossCor=T, chromosomes = args$chrom)
                                } else{
                                    exp_obj_i=ChIPQCsample(
                                        reads=file.path(samples[row,"bamReads"]),
                                        annotation = customAnnotation,
                                        runCrossCor=T, chromosomes=NULL)
                                }
                            }
                            out_path=file.path(paste0(args$out_dir,"/",samples[row,"SampleID"]))
                            runchipqc(exp_obj_i, out_path, args$facet, args$color)
                        }
                    }
                    else{
                        print(">> write sample reports given custom annotation")
                        sampleID_arr = unique(samples[args$facet])
                        for (i in 1:length(sampleID_arr[[1]])){
                            sampID=sampleID_arr[[1]][i]
                            print(paste(">>",sampID))
                            samp_design = samples[ which(samples[args$facet]==sampID), ]
                            print(samp_design)
                            if (!(is.null(args$chrom))){
                                exp_obj=ChIPQC(samp_design, 
                                    customAnnotation,
                                    chromosomes = args$chrom)
                            } else {
                                exp_obj=ChIPQC(samp_design,
                                    customAnnotation,
                                    chromosome=NULL)
                            }
                            print(">> runchipqc function")
                            out_path=file.path(paste0(args$out_dir,"/",sampID))
                            runchipqc(exp_obj, out_path, NULL, args$color)
                        }
                    }
                } else{
                    # GIVEN A GENOME WITH DEFAULT ANNOTATION

                    if (length(no_peaks) > 0){
                        samples[which(samples$Peaks %in% no_peaks),"Peaks"]=NA
                    }
                    out_path=file.path(paste0(args$out_dir))
                    
                    if((! args$indivReports) && (! args$sampleReports)){
                        print(">> write main report")
                        if (!(is.null(args$chrom))){
                            exp_obj=ChIPQC(samples, args$genome,
                                chromosomes = args$chrom)
                        } else {
                            exp_obj=ChIPQC(samples, args$genome)
                        }
                        runchipqc(exp_obj, out_path, args$facet, args$color)
                    } else if( args$indivReports){
                        print(">> write individual reports")
                        for (row in 1:nrow(samples)) {
                            print(paste(">>",samples[row]))
                            if (! (noPeakFiles)){
                                if (!(is.null(args$chrom))){
                                    exp_obj_i=ChIPQCsample(reads=file.path(samples[row,"bamReads"]),
                                        peaks = file.path(samples[row,"Peaks"]),
                                        annotation = args$genome,
                                        runCrossCor=T, chromosomes = args$chrom)
                                } else{
                                    exp_obj_i=ChIPQCsample(reads=file.path(samples[row,"bamReads"]),
                                        peaks = file.path(samples[row,"Peaks"]),
                                        annotation = args$genome,
                                        runCrossCor=T, chromosomes=NULL)
                                }
                            } else{
                                if (!(is.null(args$chrom))){
                                    exp_obj_i=ChIPQCsample(reads=file.path(samples[row,"bamReads"]),
                                        annotation = args$genome,
                                        runCrossCor=T, chromosomes = args$chrom)
                                } else{
                                    exp_obj_i=ChIPQCsample(reads=file.path(samples[row,"bamReads"]),
                                        annotation = args$genome,
                                        runCrossCor=T, chromosomes=NULL)
                                }
                            }
                            out_path=file.path(paste0(args$out_dir,"/",samples[row,"SampleID"]))
                            runchipqc(exp_obj_i, out_path, args$facet, args$color)
                        }
                    } else{
                        print(">> write sample reports")
                        sampleID_arr = unique(samples[args$facet])
                        for (i in 1:length(sampleID_arr[[1]])){
                            sampID=sampleID_arr[[1]][i]
                            print(paste(">>",sampID))
                            samp_design = samples[ which(samples[args$facet]==sampID), ]
                            if (!(is.null(args$chrom))){
                            exp_obj=ChIPQC(samp_design, args$genome,
                                chromosomes = args$chrom)
                        } else {
                            exp_obj=ChIPQC(samp_design, args$genome)
                        }
                            print(">> runchipqc function")
                            out_path=file.path(paste0(args$out_dir,"/",sampID))
                            runchipqc(exp_obj, out_path, NULL, args$color)
                        }
                    }
                }
            }
        }
    }
}
