#!/bin/bash

rna_dir="meta/lfy_rna_diffExp"
gene_bed="meta/ArabidopsisGenome/Araport11_GFF3_nonHypothetical_proteinCoding_miRNA_genes.201606.bed"
cat ${rna_dir}/Jin2020_RNASeq_DexVMock_rootCallus_t1hr_adjp0.001.txt \
	${rna_dir}/Jin2020_RNASeq_DexVMock_rootCallus_t24hr_adjp0.001.txt \
	${rna_dir}/Jin2020_RNASeq_DexVMock_rootCallus_t6hr_adjp0.001.txt \
	${rna_dir}/Winter2011_RNASeq_DexVMock_seedlings_t4hr.txt | \
	cut -f1 | grep -v "gene_id" | sort -u \
	> ${rna_dir}/lfy_dependent_genes.list

grep -f ${rna_dir}/lfy_dependent_genes.list \
    ${gene_bed} > ${rna_dir}/lfy_dependent_genes.bed
