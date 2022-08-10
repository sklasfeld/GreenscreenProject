# I. Introduction

This repository contains scripts and files used to analyze ChIP-Seq experiments for Klasfeld et al 2022. 

To generate a demo greenscreen at path `demo/data/macs2_out/inputControls/qval10/merge500bp_20inputs.txt` run the following command in the current working directory:

```
python3 scripts/greenscreenPipeline.py \
    -c Chr2 --genome_prefix TAIR10_Chr2 \
    --trimmomatic_params "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" \
    --bowtie2_index_dir demo/meta/ArabidopsisGenome/bowtie2_genome_dir \
    --chipqc_samplesheet_dir  demo/meta \
    demo/meta/input_meta.csv \
    demo/meta/ArabidopsisGenome/TAIR10_Chr2.fasta \
    demo/meta/ArabidopsisGenome/Araport11_GFF3_genes_transposons.UPDATED.201606.Chr2.gff \
    demo/meta/ArabidopsisGenome/TAIR10_Chr2_count.txt \
    demo/data
```

A docker image to run this repository can be pulled using the following command:

```
docker pull sklasfeld/greenscreen:latest
```

**A more detailed tutorial to run the greenscreen pipeline is found in TUTORIAL.pdf**

# II. Contents of This Repository

* `data` directory includes
   * *Arabidopsis* blacklist
   * *Arabidopsis* greenscreen
   * other analysis results
* `notebook` directory: jupyter notebooks used to generate some of the plots for the paper
* `meta` directory: text files used for providing context for running specific scripts 
	* See TUTORIAL.pdf and meta/README.md for more information 
	about specific files
* `scripts` directory: custom scripts used for greenscreen and ChIP-Seq analysis 
	* See TUTORIAL.pdf and scripts/README.md for more information 
	about specific files
* `demo` directory: input files to run demo script `scripts/greenscreenPipeline.py`

