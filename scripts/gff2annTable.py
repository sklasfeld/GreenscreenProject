#!/usr/bin/env python3

# -*- coding: iso-8859-15 -*-
# 2019, Samantha Klasfeld, the Wagner Lab
# the Perelman School of Medicine, the University of Pennsylvania
# Samantha Klasfeld, 10-4-2032

import os
import argparse
from argparse import RawTextHelpFormatter
import sys
import pandas as pd # used to handle dataframes/tables
import numpy as np # scientific computing python library
# libraries to translate the text
from urllib.parse import unquote
from w3lib.html import replace_entities

feature_keys=["ID", "locus_type", "Name", "Alias",
        "Dbxref", "full_name", "symbol",
        "description", "Note", "curator_summary",
        "computational_description",
        "nochangenat-description"]


def createFeatureDict(attributre_str, fkeys):
    feature_arr = attributre_str.split(";")
    feature_dict={}
    last_key=""
    for item in feature_arr:
        item_arr = item.split("=")
        if item_arr[0] in fkeys:
            feature_dict[item_arr[0]] = (
                replace_entities(
                    unquote(
                        "=".join(item_arr[1:]))))
            last_key=item_arr[0]
        else:
            feature_dict[last_key] = (
                ("%s;%s") % (
                    feature_dict[last_key], 
                    replace_entities(
                        unquote(
                            item))))
    return(feature_dict)


parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
    description="Parse attribute \
    column in GFF formatted file to create a text file for \
    gene annotations")
parser.add_argument('gff_file', help='a gff format file')
parser.add_argument('out_file', 
    help='a tab-delimited table with the following fields:\n\n' \
    "1. gene_id -Indicates the unique identifier of the gene. IDs must be unique within the scope of the GFF file. \n" \
    "2. chrom - name of the chromosome or scaffold which contains the gene\n" \
    "3. start - Start position of the gene, with sequence numbering starting at 0.\n" \
    "4. stop - End position of the feature, with sequence numbering starting at 1.\n" \
    "5. score - A floating point value.\n" \
    "6. strand - defined as + (forward) or - (reverse)\n" \
    "7. locus_type - Gene type\n" \
    "8. Name - Display name for the gene. This is the name to be displayed to the user. Unlike IDs, there is no requirement that the Name be unique within the file. \n" \
    "9. Alias - A secondary name for the gene. It is suggested that this tag be used whenever a secondary identifier for the feature is needed, such as locus names and accession numbers. Unlike ID, there is no requirement that Alias be unique within the file. \n" \
    "10. Dbxref - database cross references\n" \
    "11. full_name\n" \
    "12. symbol\n" \
    "13. description\n" \
    "14. Note\n" \
    "15. curator_summary\n" \
    "16. computational_description\n" \
    "17. nochangenat-description\n")
parser.add_argument('-f','--features', nargs="*", required=False,
    help = 'limit the output to only contain specific features from gff file')
args = parser.parse_args()

# import gff table
gff_df = pd.read_csv(args.gff_file, sep="\t",
    comment='#',
    names=["chrom","source","feature","start",
    "end","score","strand","frame","attribute"])

# filture for gene features
if args.features:
    for feat in args.features:
        gff_df = gff_df.loc[gff_df["feature"]==feat,:].copy()

# generate the table
out_table = open(args.out_file, 'w')
out_table.write("gene_id\tchrom\tstart\tstop\tscore\tstrand\t%s\n" %
    ("\t".join(feature_keys[1:])))
for index, row in gff_df.iterrows():
    gene_location=("%s\t%s\t%s\t%s\t%s" % 
        (row["chrom"], row["start"]-1, 
            row["end"], row["score"],
            row["strand"]))
    feature_dict = createFeatureDict(row["attribute"],
        feature_keys)
    out_table.write("%s\t%s\t%s\t%s\t%s\t%s" %
        (feature_dict["ID"],
            row["chrom"], row["start"]-1, 
            row["end"], row["score"],
            row["strand"]))
    for possible_key in feature_keys[1:]:
        if possible_key in feature_dict:
            out_table.write("\t%s" % (feature_dict[possible_key]))
        else:
            out_table.write("\t")
    out_table.write("\n")