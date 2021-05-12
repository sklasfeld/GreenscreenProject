#!/usr/bin/env python3

# -*- coding: iso-8859-15 -*-
# 2021, Samantha Klasfeld, the Wagner Lab
# the Perelman School of Medicine, the University of Pennsylvania
# Samantha Klasfeld, 5-11-2021

import os
import argparse
import sys
import subprocess
from collections import defaultdict
import pandas as pd
import numpy as np
import pyBigWig

parser = argparse.ArgumentParser(description="Get coverage matrix \
        where rows are each region and columns are each sample.")
parser.add_argument('infile', help='comma-delimited table that \
    contains columns: `sample_name` and `mapping_file`. The \
    `sample_name` contains the ID for each respective `mapping_file`. \
    The `mapping_file` should be in bigwig format')
parser.add_argument('peaks_bed', help='bedfile of all the regions we \
    want coverage for')
parser.add_argument('-o','--outdir_path', help='path to output files', \
    default = "")
parser.add_argument('-m','--out_matrix', help='file name for matrix to output', \
    default = "coverage_matrix.csv")
parser.add_argument('-r', '--round2ints', help='counts must be in integers \
    (Default: FALSE)', action='store_true')
parser.add_argument('-v', '--verbose', help='if set, print log to standard out',
 action='store_true')
args = parser.parse_args()

# Get Full Paths
outdir = ""
if len (args.outdir_path) > 0 :
    if not os.path.isdir(args.outdir_path):
        mkdir_out = ("mkdir %s" % args.outdir_path)
        os.system(mkdir_out)
    outdir = ("%s/" % os.path.abspath(args.outdir_path))


# Get Input information
if not os.path.isfile(args.infile):
    sys.exit("ERROR: cannot find input file: %s" % args.infile)

meta_data = pd.read_csv(args.infile, sep=",")
meta_data_cols = ["sample_name","mapping_file"]
for md_col in meta_data_cols:
    if not md_col in list(meta_data.columns):
        sys.exit(("ERROR: infile %s must contain column %s") % 
            (args.infile, md_col))
for bamF in list(meta_data["mapping_file"].unique()):
    if not os.path.exists(bamF):
        sys.exit(("ERROR: BAM file %s does not exist.") 
            % (bamF))


# get bed dataframe
bed_cols_list = ["chrom", "start", "stop"]
bed_df = pd.read_csv(args.peaks_bed, sep="\t", header=None, \
        names=bed_cols_list, index_col=False, usecols=range(0,3),
        dtype={"chrom":str, \
        "start":np.int64, "stop":np.int64})

final_matrix_df = bed_df.copy()

# use mapping files to get coverage and matrix

read_counts=[]

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
                sys.stderr.write(("Chromosome %s is not in the bigwig file.") % 
                    (bed_row["chrom"]))
            else:
                sys.stderr.write(("Region %i-%i goes beyond Chromosome %s of length %i\n") % 
                    (int(bed_row["start"]), int(bed_row["stop"]),
                       bed_row["chrom"],
                        bw.chroms(bed_row["chrom"])))
            merge_df.loc[bed_index, "num_feature_overlap"] = 0
            merge_df.loc[bed_index, "length"] = bed_row["stop"] - bed_row["start"]

    
    final_matrix_df[baseOfile] = merge_df["num_feature_overlap"]

if args.round2ints:
    final_matrix_df[samples] = np.round(final_matrix_df[samples]).astype(int)

final_matrix_file = ("%s%s" % (outdir,args.out_matrix))
final_matrix_df.to_csv(final_matrix_file, sep=",", index=False)
