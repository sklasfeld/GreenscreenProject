#!/bin/python3

# -*- coding: iso-8859-15 -*-
# 2021, Samantha Klasfeld, the Wagner Lab
# the Perelman School of Medicine, the University of Pennsylvania
# Samantha Klasfeld, 9-7-2021

import sys
import argparse
from urllib.parse import unquote

parser = argparse.ArgumentParser(description="Translate \
    a file that contains unicode")
parser.add_argument('infile', help='in file containing percent unicode')
parser.add_argument('outfile', help='out file without unicode')
args = parser.parse_args()

# Write to file
new_gff = open(args.outfile, 'w')

# Opening file
org_gff = open(args.infile, 'r')

with open(args.infile) as org_gff:
    Lines = org_gff.readlines()
    for lines in Lines:
        new_gff.write("%s\n" % unquote(lines).strip())

# Closing files
org_gff.close()
new_gff.close()
