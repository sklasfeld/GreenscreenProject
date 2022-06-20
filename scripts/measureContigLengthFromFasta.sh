#!/bin/bash

# set `memory_gt_ram` to 0 if more ram is available than long-term memory
# keep `memory_gt_ram` set to 1 if otherwise
memory_gt_ram=1

# from decompressed fastq get table output
contig_size () {
  awk '{ \
	if($1~">"){ \
		if(NR>1){ \
		print chrom"\t"n_nuc} \
		n_nuc=0; \
		chrom=substr($1,2) \
		}else{ \
		n_nuc=n_nuc + length($0)}} \
	END{print chrom"\t"n_nuc}' \
	$1 > $2
}


if [ "$#" -lt 2 ]; then
	echo -e "\n> Goal: Calculate genome contig sizes"
	echo -e "\n> Command Format:"
	echo -e "\nbash /path/to/measureContigLengthFromFasta.sh [genome:REQUIRED] [genomesize_table_output:REQUIRED] [memory_gt_ram=1]"
	echo -e "\n> Required Parameters:"
	echo -e "\n\t* genome: A FASTA-formatted file containing genome nucleotide sequences\n"
	echo -e "\n\t* genomesize_table_output: path to output file which table will be printed to\n"
	echo -e "\n> Optional Parameters:"
	printf "\n\t* memory_gt_ram: Set to 1 or 0 based on computational limiting factors (Default=1)."
	printf " Set this to 0 if more ram is available than long-term memory. Otherwise keep set to 1.\n\n"
elif [ "$#" -gt 2 ] && [ "$3" != 1 ] && [ "$3" != 0 ]; then
	printf "\nError: Set third parameter to either 0 or 1. "
	echo -e "Set to 0 if more ram is available than long-term memory. Otherwise keep set to 1."
else
	if [ "$#" -gt 2 ]; then
		memory_gt_ram=$3
	fi

	fasta_f=$1
	genomesize_table_output=$2
	
	if [[ ${memory_gt_ram} -eq 0 ]]; then
		if [[ $fasta_f == *.gz ]]; then
			gunzip -c ${fasta_f} | \
				awk '{ \
					if($1~">"){ \
						if(NR>1){ \
						print chrom"\t"n_nuc} \
						n_nuc=0; \
						chrom=substr($1,2) \
						}else{ \
						n_nuc=n_nuc + length($0)}} \
					END{print chrom"\t"n_nuc}' > \
				${genomesize_table_output}
		else
			contig_size ${fasta_f} ${genomesize_table_output}
		fi
	else
		rezip=0
		if [[ $fasta_f == *.gz ]]; then
			rezip=1
			# code below is best if you 
			# have more memory than RAM
			gunzip ${fasta_f}
			suffix=".gz"
			fasta_f=${fasta_f%"$suffix"}
		fi

		contig_size ${fasta_f} ${genomesize_table_output}
		if [[ $rezip -eq 1 ]]; then
			gzip ${fasta_f}
		fi
	fi
fi
