# Main Directory

* `meta/chip_controls_fragsize_nreps.csv`

	* Respective ChIP-Seq sample names, number of ChIP-Seq replicates,
	ChIP-Seq sample fragment sizes, ChIP-Seq control names, 
	number of ChIP-Seq control replicates

	* Each row contains:

		* column 1: ChIP-Seq sample name

		* column 2: Number of replicates of ChIP-Seq sample

		* column 3: ChIP-Seq sample fragment size

		* column 4: ChIP-Seq control name

		* column 5: Number of replicates of ChIP-Seq control

* `meta/chip_controls_sampleIDs.list`

	* List delimited by new lines that contains all the ChIP-Seq sample and control names

* `meta/chip_readsize_fragsize.csv`

	* Respective ChIP-Seq replicate names, read size, and fragment sizes

	* Each row contains:

		* column 1: ChIP-Seq sample replicate names

		* column 2: ChIP-Seq sample read size

		* column 3: ChIP-Seq sample fragment size

* `meta/ChIPseqSampleSheet.csv`

	* Sample sheet for ChIPQC for all ChIP-Seq experiments

	* Each row contains:

		* `SampleID`: ChIP-Seq sample name

		* `bamReads`: File path to respective mapped reads with MAPQ>=30

* `meta/chip_trueRep_bigwigs.csv`

	* Table which matches ChIP-seq replicate names (`sample_name`)
	to their respective bigwig paths (`mapping_file`)

* `meta/chip_trueReps_colorshapeLabels.csv`

	* Table which matches ChIP-Seq replicates to specific colors and shapes
	to label the heatmap exported by `scripts/readCorrelationPlot.py`

* `meta/chip_trueReps_expectedCluster.csv`

	* Table which matches ChIP-Seq replicates to expected clusters.
	This is used to calculate the rand-index between expected verses
	unsupervised hierarchical clusters in `scripts/readCorrelationPlot.py` 

* `meta/input_readsizes.csv`

	* Respective input letter and read sizes for each of the random inputs 
	used to call Artificial Signals

	* Each row contains:

		* column 1: Input sample letter (A-T)

		* column 2: Input read size


* `meta/InputseqSampleSheet.csv`

	* Sample sheet for ChIPQC for random inputs used to call Artificial Signals

	* Each row contains:

		* `SampleID`: Input sample name

		* `bamReads`: File path to respective mapped reads with MAPQ>=30

* `meta/chip_trueRep_bigwigs.csv`

	* table that contains columns: `sample_name`,`mapping_file`,
                        (and `seq_depth`). `seq_depth` is only necessary if `normalize` parameter
                        is set. To calculate sequencing depth for each sample use samtools.

# lfy_rna_diffExp

Differential expression metrics of genes that displayed rapid changes in gene 
expression after LFY-GR activation

* `meta/lfy_rna_diffExp/Jin2020_RNASeq_DexVMock_rootCallus_t1hr_adjp0.001.txt`

	* Source: Jin 2021

	* Sample: Differential expression was measured between root explants 
	treated with 5 µM dexamethasone for 1 h before the end of the 5-day 
	incubation on CIM plates to those treated with 5 µM mock solution for
	1 h before the end of the 5-day incubation on CIM plates

	* DESeq2 q-value <= 0.01

* `meta/lfy_rna_diffExp/Jin2020_RNASeq_DexVMock_rootCallus_t6hr_adjp0.001.txt`

	* Source: Jin 2021

	* Sample: Differential expression was measured between root explants 
	treated with 5 µM dexamethasone for 6 h before the end of the 5-day 
	incubation on CIM plates to those treated with 5 µM mock solution for
	6 h before the end of the 5-day incubation on CIM plates

	* DESeq2 q-value <= 0.01
	
* `meta/lfy_rna_diffExp/Jin2020_RNASeq_DexVMock_rootCallus_t24hr_adjp0.001.txt`

	* Source: Jin 2021

	* Sample: Differential expression was measured between root explants 
	treated with 5 µM dexamethasone for 24 h before the end of the 5-day 
	incubation on CIM plates to those treated with 5 µM mock solution for
	24 h before the end of the 5-day incubation on CIM plates

	* DESeq2 (q-value <= 0.01)

* `meta/lfy_rna_diffExp/Winter2011_RNASeq_DexVMock_seedlings_t4hr.txt`

	* Source: Winter 2011

	* Sample: Differential expression was measured between shoot apices 
	from 9-day-old 35S:LFY-GR seedlings treated with 10 μM dexamethasone
	in 0.1% ethanol, for 4 hr as described (Wagner 1999) to 
	9-day old *Ler* seedlings given the same treatment

	* LIMMA (FDR<0.05 and |fold-change|>1.5)
