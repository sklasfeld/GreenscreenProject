#!/usr/bin/env python3

# -*- coding: iso-8859-15 -*-
# 2019, Samantha Klasfeld, the Wagner Lab
# the Perelman School of Medicine, the University of Pennsylvania
# Samantha Klasfeld, 12-9-2019

import os
import argparse
import sys
import subprocess
from collections import defaultdict
import pandas as pd
import numpy as np
import pyBigWig

parser = argparse.ArgumentParser(
    description="Given specific regions and multiple bigwig files, \
    generate a signal matrix where each row represents a specific \
    region and each column represents a sequencing experiment")
parser.add_argument('infile', help='comma-delimited table that \
    matches sample names to bigwig files. The sample \
    name column should be labeled `sample_name`, and column \
    containing the bigwig files should be labeled `mapping_file`. \
    Optional: Add a column `seq_depth` containing the total read \
    depth of a sample if `--normalize` parameter is set to \
    `CPM`or `RPKM`.')
parser.add_argument('peaks_bed', help='BED-formatted file \
    which contains the regions of interest')
parser.add_argument('-o','--outdir_path', help='path to output files', \
    default = "")
parser.add_argument('-bp','--bedtools_path', help='path to bedtools', \
    default = "")
parser.add_argument('-m','--out_matrix', help='file name for matrix to output', \
    default = "coverage_matrix.csv")
parser.add_argument('-n','--normalize',
    help='Normalize by CPM, RPKM, or scaling. Note that scaling means \
    down-scaling to the sample with the lowest amount of reads.')
parser.add_argument( '-s', '--scale', help='scaling factor (eg. to get RPKM, \
    use 1e6). Only set when `--normalize` flag is set. Default: 1e6', \
    type=float, default=1e6)
parser.add_argument('-r', '--round2ints', help='counts must be in integers \
    (Default: FALSE)', action='store_true')
parser.add_argument('-v', '--verbose', help='if set, print log to standard out',
 action='store_true')
args = parser.parse_args()

if args.normalize:
    args.normalize = args.normalize.upper()
    if (args.normalize != "RPKM" and
        args.normalize != "CPM" and
        args.normalize != "SCALING"):
        sys.exit("Can only normalize by RPKM, CPM, or scaling")

# Local Methods
def wc(file_name, verbose=False):
    if (not os.path.isfile(file_name)):
        if verbose:
            sys.stdout.write("Count: 0\n")
            sys.stdout.flush()
        return 0
    else:
        cmd = ("wc -l %s" % file_name)
        if os.path.splitext(file_name)[1] == ".gz":
            cmd = ("gunzip -c %s | wc -l" % file_name)
        if verbose:
            sys.stdout.write("%s\n" % cmd)
            sys.stdout.flush()
        ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
        output = ps.communicate()[0]
        break_out = output.split()
        if verbose:
            sys.stdout.write("Count: %i\n" % int(break_out[0]))
            sys.stdout.flush()
        return int(break_out[0])

# Get Full Paths
outdir = ""
if len (args.outdir_path) > 0 :
    if not os.path.isdir(args.outdir_path):
        mkdir_out = ("mkdir %s" % args.outdir_path)
        os.system(mkdir_out)
    outdir = ("%s/" % os.path.abspath(args.outdir_path))

bedtools_p = ""
if len (args.bedtools_path) > 0 :
    bedtools_p = ("%s/" % os.path.abspath(args.bedtools_path))

# Get Input information
if not os.path.isfile(args.infile):
    sys.exit("ERROR: cannot find input file: %s" % args.infile)

meta_data = pd.read_csv(args.infile, sep=",")
meta_data_cols = ["sample_name","mapping_file"]
if args.normalize:
    if (args.normalize == "RPKM" or
        args.normalize == "CPM"):
        meta_data_cols.append("seq_depth")
for md_col in meta_data_cols:
    if not md_col in list(meta_data.columns):
        sys.exit(("ERROR: infile %s must contain column %s") % 
            (args.infile, md_col))
for mapF in list(meta_data["mapping_file"].unique()):
    if not os.path.exists(mapF):
        sys.exit(("ERROR: mapping file %s does not exist.") 
            % (mapF))
    # check that files are bigwig format    
    else:
        file_extension = os.path.splitext(mapF)[1]
        if (not (file_extension == ".bigwig" or \
        file_extension == ".bw")):
            sys.exit(("ERROR: The mapping file (%s) " +
            "does not seem to be in bigwig format." +
            "The suffix is neither `.bw` nor `.bigwig`") 
            % (mapF))

# get bed dataframe
bed_cols_list = ["chrom", "start", "stop"]
bed_df = pd.read_csv(args.peaks_bed, sep="\t", header=None, \
        names=bed_cols_list, index_col=False, usecols=range(0,3),
        dtype={"chrom":str, \
        "start":np.int64, "stop":np.int64})

final_matrix_df = bed_df.copy()

# use mapping files to get coverage and matrix

read_counts=[]
if args.normalize:
    read_counts=list(meta_data["seq_depth"])

for index, row in meta_data.iterrows():
    
    bfile = row["mapping_file"]
    baseOfile = row["sample_name"]
    if args.verbose:
        sys.stdout.write("Analyzing %s (%s)\n" % (bfile,baseOfile))
        sys.stdout.flush()
    
    bw = pyBigWig.open(bfile)
    if not bw.isBigWig(): # check bigwig format
        sys.exit("ERROR: %s is not a bigwig file!" % bfile)
    merge_df = bed_df.copy()
    for bed_index, bed_row in merge_df.iterrows():
        if (bed_row["chrom"] in bw.chroms().keys() and 
            int(bed_row["start"]) < int(bw.chroms(bed_row["chrom"]))):
            if int(bed_row["stop"]) >= int(bw.chroms(bed_row["chrom"])):
                sys.stderr.write(("Region %i-%i overlaps the end of " +
                    "Chromosome %s of length %i\n") % 
                (int(bed_row["start"]), int(bed_row["stop"]),
                    bed_row["chrom"],
                    bw.chroms(bed_row["chrom"])))
                bed_row["stop"] = bw.chroms(bed_row["chrom"])
            merge_df.loc[bed_index, "num_feature_overlap"] = \
                bw.stats(bed_row["chrom"],
                bed_row["start"],
                bed_row["stop"],
                type = "max")[0]
            merge_df.loc[bed_index, "length"] = bed_row["stop"] - bed_row["start"]
        else:
            sys.stderr.write(("%s\n") % (bed_row["start"]))
            sys.stderr.write(("%s\n") % (bw.chroms(bed_row["chrom"])))
            if (bed_row["chrom"] not in bw.chroms().keys()):
                if args.verbose:
                    sys.stderr.write(("Chromosome %s is not in the bigwig file.") % 
    	                (bed_row["chrom"]))
            else:
                sys.stderr.write(("Region %i-%i goes beyond Chromosome %s of length %i\n") % 
                    (int(bed_row["start"]), int(bed_row["stop"]),
                       bed_row["chrom"],
                        bw.chroms(bed_row["chrom"])))
            merge_df.loc[bed_index, "num_feature_overlap"] = 0
            merge_df.loc[bed_index, "length"] = bed_row["stop"] - bed_row["start"]
    if args.normalize == "RPKM": 
        final_matrix_df[baseOfile] = merge_df["num_feature_overlap"] / \
            (merge_df["length"] * scaling_factor)
    elif args.normalize == "CPM":
        final_matrix_df[baseOfile] = merge_df["num_feature_overlap"] / \
            (scaling_factor)
    else:
        final_matrix_df[baseOfile] = merge_df["num_feature_overlap"]
if args.normalize == "SCALING":
    min_rc = float(np.min(read_counts))
    for index, row in meta_data.iterrows():
        sample_name = row["sample_name"]
        read_count = row["seq_depth"]
        scale_factor=read_count/min_rc
        final_matrix_df[sample_name] = \
            final_matrix_df[sample_name] / scale_factor

if args.round2ints:
    final_matrix_df[samples] = np.round(final_matrix_df[samples]).astype(int)

final_matrix_file = ("%s%s" % (outdir,args.out_matrix))
final_matrix_df.to_csv(final_matrix_file, sep=",", index=False)
