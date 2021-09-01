#!/bin/bash

awk -F"\t" \
    '$3=="gene"{ \
        locus_type=0; \
        hypoth=0; \
        nmeta=split($9,meta,";"); \
        for(i=1;i<=nmeta;i++){ \
            if(meta[i]~"ID="){ \
                split(meta[i],id_arr,"="); \
                id=id_arr[2] \
            }if(meta[i]~"locus_type=protein_coding" || \
                meta[i]~"locus_type=mirna"){ \
                    locus_type=1 \
            }if(meta[i]~"Note=" && meta[i]~"hypothetical protein"){ \
                hypoth=1 \
            } \
        } if(locus_type==1 && hypoth!=1){ \
            print $1"\t"$4-1"\t"$5"\t"id"\t"$6"\t"$7}}' \
    meta/ArabidopsisGenome/Araport11_GFF3_genes_transposons.201606.gff \
    > meta/ArabidopsisGenome/Araport11_GFF3_nonHypothetical_proteinCoding_miRNA_genes.201606.bed
