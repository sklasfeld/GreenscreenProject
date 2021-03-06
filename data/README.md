# Contents of this Directory
This directory contains tables and files used to generate some of the plots for the paper
* arabidopsis_blacklist_20inputs.bed
  * blacklist regions were generated using the [blacklist tool](https://github.com/Boyle-Lab/Blacklist) (Amemiya 2019)
  * the blacklist tool was manually modified  [here](https://github.com/Boyle-Lab/Blacklist/blob/master/blacklist.cpp#L469)
    * original line of code: `if(miss < (binSize/binOverlap + 200)) { // bridge over adjacent bins plus 100 * 200 = 20kb`
    * modified line of code: `if(miss < (binSize/binOverlap + 50)) { // bridge over adjacent bins plus 100 * 200 = 20kb`
  * Input to Blacklist Tool:
    * output files from UMAP (v1.1.1)
      * Note that the repo for this tool was originally located at the following (link)[https://bitbucket.org/hoffmanlab/proj/bismap]
      * However, the repo was moved to this (link)[https://github.com/hoffmangroup/umap]
    * Twenty input-seq samples from *Arabidopsis thaliana* were cleaned and mapped (MAPQ>3=0). See Table S2 in the Greenscreen Paper.
* arabidopsis_greenscreen_20inputs.bed
  * Green screen regions were generated using the green screen pipeline.
    1. MACS2 v2.2.7.1 (—keepdup “auto” —nomodel –extsize [read length] —broad --nolambda  -g 101274395) was run on each of the twenty inputs (Table S2 in the paper). 
    2. Peaks with average q-values (column 9 in the broadPeak output file) above 10-10 were filtered out. 
    3. Peaks called from the individual inputs were concatenated
    4. All peaks within a 5kb distance were merged. 
    5. Regions that did not overlap at least half of the inputs analyzed were removed
* fd_pooled_chipComparisons.tsv
  * This table is in reference to ChIP-seq peaks called by MACS2 v2.2.7.1 (--keepdup "auto" --nomodel --extsize 250  -g 101274395) from FD replicates (Zhu *et al* 2020) that were down-sampled and pooled (FD_W)
  * FD_W is compared to:
    * green screen peaks (data/arabidopsis_greenscreen_20inputs.bed)
    * ChIP-seq peaks called by MACS2 v2.2.7.1 (--keepdup "auto" --nomodel --extsize 220  -g 101274395) from FD replicates (Zhu *et al* 2020) that were down-sampled and pooled (FD_S) using a mock sample generated from the same experiment as the control.
  * different parameters were tested on FD_W to see how it changes the peak overlap
* lfy_pooled_chipComparisons.tsv
  * This table is in reference to ChIP-seq peaks called by MACS2 v2.2.7.1 (--keepdup "auto" --nomodel --extsize 150  -g 101274395) from LFY replicates (Jin *et al* 2021) that were down-sampled and pooled (LFY_W)
  * LFY_W is compared to:
    * green screen peaks (data/arabidopsis_greenscreen_20inputs.bed)
    * LFY ChIP-chip peaks from Winter *et al* 2011 
  * different parameters were tested on LFY_W to see how it changes the peak overlap
* callArtifSignalMetrics.tsv
  * Metrics table to compare how different Arabidopsis artificial signals (AS) masks alter downstream analysis results
  * Features:
    * AS_mask_type: AS Mask Type
      * None: No mask
      * blacklist: reads removed that overlapped with blacklist regions
      * greenscreen: peaks removed that overlapped with green screen regions
    * average_MACS2_qvalue: Green screen parameter to call peaks in inputs that have at least this set mean -log10 q-value of a peak region.
      * For example, if this parameter is set to 10 then peaks called by MACS2 are filtered to have a mean q-value <= 10^-10.
    * merge_distance: Blacklist and greenscreen parameter to merge two significant regions within this maximum distance.
      * For example, if this papameter is set to 1000 then significant regions called by either the blacklist tool or MACS2 for green screen that are within 1kb apart will become a single region. 
    * num_inputs_used: Number of inputs imported into either the blacklist tool or MACS2 for green screen
    * num_inputs_available: To test the number of random inputs needs for identifying artificial signals, subsamples were normally chosen from the 20 random inputs that were analyzed. However, there are only the four inputs to match the ChIP experiments from Zhu et al 2021 in the downstream analysis. To compare the results of the downstream analysis between random and related input samples, we performed analysis using 3 random/related input samples chosen from 4 random/related input samples.
    * min_samples_called: Green screen parameter sets the minimum number of input samples' peaks must overlap for a region to be labeled as part of the green screen
      * regions where less than `min_samples_called` samples have peaks that overlap will not be included in the green screen
    * inputs_used: The names of the inputs used. 
      * information about inputs A-T can be found in Table S2.
      * TFL1_Input_R1, TFL1_Input_R2, TFL1_Input_R3, TFL1_Input_R4 
	are the inputs used in et al Zhu 2021.
    * genome_cov: number of basepairs that the AS regions mask out
    * 100per_genes_coverage: the number of genes 100% masked out by the AS mask
    * perc_lfy_chipseq_in_chipchip: percent of LFY callus peaks given MASK overlap with LFY chip-chip seedlings
    * c(Yk=2, Y’k=2): rand-index values to measure k=2 clustering hypothesis
    * c(Yk=3, Y’k=3): rand-index values to measure k=3 clustering hypothesis
    * c(Zk=4, Z’k=4): rand-index values to measure k=4 clustering hypothesis
* ChIPQCreport
  * directory containing the results of ChIPQC on inputs and ChIP-seq experiments used in this paper
  * see README.md in this directory for more information
