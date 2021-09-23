#!/bin/bash

if [ "$#" -ne 1 ]; then
	echo -e "\n> Goal: Calculate genome contig sizes"
	echo -e "\n> Command Format:"
	echo -e "\nbash /path/to/measureContigLengthFromFasta.sh [genome]"
	echo -e "\n> Parameters:"
	echo -e "\n\t* genome: A FASTA-formatted file containing genome nucleotide sequences\n"
else

	fasta_f=$1

	if [[ $fasta_f == *.gz ]]; then
		gunzip ${fasta_f}
		suffix=".gz"
		fasta_f=${fasta_f%"$suffix"}
	fi

	awk '{ \
    	if($1~">"){ \
        	if(NR>1){ \
            	print chrom"\t"n_nuc} \
        	n_nuc=0; \
        	chrom=substr($1,2) \
        	}else{ \
            	n_nuc=n_nuc + length($0)}} \
	    END{print chrom"\t"n_nuc}' \
	    ${fasta_f}
fi
