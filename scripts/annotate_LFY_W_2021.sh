#!/bin/bash

# paths to essential directories
summits_dir="data/macs2_out/chipPeaks/gsMask_qval10"
out_dir="data/annotations/LFY_Jin_2021"
rna_dir="meta/lfy_rna_diffExp"

python3 scripts/ChIP_Annotation LFY_W_2021 \
    ${out_dir}/ \
    ${summits_dir}/LFY_W_2021_summits.bed \
    meta/ArabidopsisGenome/Araport11_GFF3_true_protein_miRNA_genes.bed \
    -n ${summits_dir}/LFY_W_2021_summits.narrowPeak \
    -tss 4000 -tts 0 \
    -gt meta/ArabidopsisGenome/Araport11_GFF3_genesAttributes.tsv \
    ${rna_dir}/Winter2011_RNASeq_DexVMock_seedlings_t4hr.txt \
    ${rna_dir}/Jin2020_RNASeq_DexVMock_rootCallus_t1hr_adjp0.001.txt \
    ${rna_dir}/Jin2020_RNASeq_DexVMock_rootCallus_t24hr_adjp0.001.txt \
    ${rna_dir}/Jin2020_RNASeq_DexVMock_rootCallus_t6hr_adjp0.001.txt \
    -gtn "geneMeta" "RNA_DE_seed" \
    "RNA2_DE_callus_t1hr" "RNA2_DE_callus_t6hr" "RNA2_DE_callus_t24hr" \
    -gtf "geneMeta:Name,Alias,Note,locus_type,description" \
    "RNA_DE_seed:log2FoldChange,padj" "RNA2_DE_callus_t1hr:log2FoldChange,padj" \
    "RNA2_DE_callus_t6hr:log2FoldChange,padj" "RNA2_DE_callus_t24hr:log2FoldChange,padj" \
    -r2 ${rna_dir}/Winter2011_RNASeq_DexVMock_seedlings_t4hr.txt \
    ${rna_dir}/Jin2020_RNASeq_DexVMock_rootCallus_t1hr_adjp0.001.txt \
    ${rna_dir}/Jin2020_RNASeq_DexVMock_rootCallus_t24hr_adjp0.001.txt \
    ${rna_dir}/Jin2020_RNASeq_DexVMock_rootCallus_t6hr_adjp0.001.txt
