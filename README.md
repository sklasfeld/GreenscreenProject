# I. Contents of This Repository
## A. `Data` Directory
The `Data` directory contains tables and files used to generate some of the plots for the paper
* data/arabidopsis_blacklist_20inputs.bed
 * blacklist regions were generated using the [blacklist tool](https://github.com/Boyle-Lab/Blacklist) (Amemiya 2019)
 * the blacklist tool was manually modified  [here](https://github.com/Boyle-Lab/Blacklist/blob/master/blacklist.cpp#L469)
   * original line of code: `if(miss < (binSize/binOverlap + 200)) { // bridge over adjacent bins plus 100 * 200 = 20kb`
   * modified line of code: `if(miss < (binSize/binOverlap + 50)) { // bridge over adjacent bins plus 100 * 200 = 20kb`
 * Input to Blacklist Tool:
   * output files from UMAP (v1.1.1)
    * Note that the repo for this tool was originally located at the following (link)[https://bitbucket.org/hoffmanlab/proj/bismap]
    * However, the repo was moved to this (link)[https://github.com/hoffmangroup/umap]
  * Twenty input-seq samples from *Arabidopsis thaliana* were cleaned and mapped (MAPQ>3=0). See Table S2 in the Greenscreen Paper.
* data/arabidopsis_greenscreen_20inputs.bed
 * Green screen regions were generated using the green screen pipeline.
   1. MACS2 v2.2.7.1 (—keepdup “auto” —nomodel –extsize [read length] —broad --nolambda  -g 101274395) was run on each of the twenty inputs (Table S2 in the paper). 
    2. Peaks with average q-values (column 9 in the broadPeak output file) above 10-10 were filtered out. 
    3. Peaks called from the individual inputs were concatenated
    4. All peaks within a 5kb distance were merged. 
    5. Regions that did not overlap at least half of the inputs analyzed were removed
 * data/fd_pooled_chipComparisons.tsv
   * This table is in reference to ChIP-seq peaks called by MACS2 v2.2.7.1 (--keepdup "auto" --nomodel --extsize 250  -g 101274395) from FD replicates (Zhu *et al* 2020) that were down-sampled and pooled (FD_W)
   * FD_W is compared to:
    * green screen peaks (data/arabidopsis_greenscreen_20inputs.bed)
    * ChIP-seq peaks called by MACS2 v2.2.7.1 (--keepdup "auto" --nomodel --extsize 220  -g 101274395) from FD replicates (Zhu *et al* 2020) that were down-sampled and pooled (FD_S) using a mock sample generated from the same experiment as the control.
   * different parameters were tested on FD_W to see how it changes the peak overlap
 * data/lfy_pooled_chipComparisons.tsv
   * This table is in reference to ChIP-seq peaks called by MACS2 v2.2.7.1 (--keepdup "auto" --nomodel --extsize 150  -g 101274395) from LFY replicates (Jin *et al* 2021) that were down-sampled and pooled (LFY_W)
   * LFY_W is compared to:
    * green screen peaks (data/arabidopsis_greenscreen_20inputs.bed)
    * LFY ChIP-chip peaks from Winter *et al* 2011 
   * different parameters were tested on LFY_W to see how it changes the peak overlap
* data/ChIPQCreport/ : directory containing the results of ChIPQC on inputs and ChIP-seq experiments used in this paper
  * 20inputs_blMask
   * ChIPQC Results for 20 random inputs (see Table S2 in the paper) after using samtools to remove reads that overlap blacklist (data/arabidopsis_blacklist_20inputs.bed)
  * 20inputs_gsMask
   * ChIPQC Results for 20 random inputs (see Table S2 in the paper) after using samtools to remove reads that overlap greenscreen (data/arabidopsis_greenscreen_20inputs.bed)
  * 20inputs_noMask
   * ChIPQC Results for 20 random inputs (see Table S2 in the paper) 
  * chip_blMask_noDups_lfyWInputControl_extFragSize_qval10
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads (without duplicates) after blacklist (data/arabidopsis_blacklist_20inputs.bed) masking
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Jin 2021 were down-sampled and pooled to use as the MACS2 control
   * Peaks with summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included
  * chip_blMask_noDups_tfl1WInputControl_extFragSize_qval10:
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads (without duplicates) after blacklist (data/arabidopsis_blacklist_20inputs.bed) masking
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Zhu 2020 were down-sampled and pooled to use as the MACS2 control
   * Peaks with summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included
  * chip_blMask_wiDups_lfyWInputControl_extFragSize_qval10: 
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads after blacklist (data/arabidopsis_blacklist_20inputs.bed) masking
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Jin 2021 were down-sampled and pooled to use as the MACS2 control
   * Peaks with summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included
  * chip_blMask_wiDups_noControl_extFragSize_qval10: 
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads (without duplicates) after blacklist (data/arabidopsis_blacklist_20inputs.bed) masking
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Nothin was set as the MACS2 control
   * Peaks with summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included
  * chip_blMask_wiDups_respectiveMockControls_extFragSize_qval10:
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads (without duplicates) after blacklist (data/arabidopsis_blacklist_20inputs.bed) masking
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Mock samples were down-sampled and pooled to use as the MACS2 control for each ChIP respectively (see Table S3 in the paper)
   * Peaks with summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included
  * chip_blMask_wiDups_tfl1WInputControl_extFragSize_qval10: 
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads (without duplicates) after blacklist (data/arabidopsis_blacklist_20inputs.bed) masking
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Zhu 2020 were down-sampled and pooled to use as the MACS2 control
   * Peaks with summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included
  * chip_gsMask_noDups_lfyWInputControl_extFragSize_qval10:
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads (without duplicates)
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Jin 2021 were down-sampled and pooled to use as the MACS2 control
   * Peaks that overlap green screen regions (data/arabidopsis_greenscreen_20inputs.bed) OR had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included
  * chip_gsMask_noDups_tfl1WInputControl_extFragSize_qval10: 
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads (without duplicates)
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Zhu 2020 were down-sampled and pooled to use as the MACS2 control
   * Peaks that overlap green screen regions (data/arabidopsis_greenscreen_20inputs.bed) OR had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included
  * chip_gsMask_wiDups_lfyWInputControl_extFragSize_qval10:
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Jin 2021 were down-sampled and pooled to use as the MACS2 control
   * Peaks that overlap green screen regions (data/arabidopsis_greenscreen_20inputs.bed) OR had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
  * chip_gsMask_wiDups_noControl_extFragSize_qval10:
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Nothing was set as the MACS2 control
   * Peaks that overlap green screen regions (data/arabidopsis_greenscreen_20inputs.bed) OR had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
  * chip_gsMask_wiDups_respectiveMockControls_extFragSize_qval10: 
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Mock samples were down-sampled and pooled to use as the MACS2 control for each ChIP respectively (see Table S3 in the paper)
   * Peaks that overlap green screen regions (data/arabidopsis_greenscreen_20inputs.bed) OR had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
  * chip_gsMask_wiDups_tfl1WInputControl_extFragSize_qval10:
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Zhu 2020 were down-sampled and pooled to use as the MACS2 control
   * Peaks that overlap green screen regions (data/arabidopsis_greenscreen_20inputs.bed) OR had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
  * chip_noMask_noDups_lfyWInputControl_extFragSize_qval10:
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads (without duplicates)
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Jin 2021 were down-sampled and pooled to use as the MACS2 control
   * Peaks that had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
  * chip_noMask_noDups_tfl1WInputControl_extFragSize_qval10:
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads (without duplicates)
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Zhu 2020 were down-sampled and pooled to use as the MACS2 control
   * Peaks that had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
  * chip_noMask_wiDups_lfyWInputControl_extFragSize_qval10:
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Jin 2021 were down-sampled and pooled to use as the MACS2 control
   * Peaks that had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
  * chip_noMask_wiDups_noControl_extFragSize_qval10: 
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Nothing was used as the MACS2 control
   * Peaks that had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
  * chip_noMask_wiDups_respectiveMockControls_extFragSize_qval10:
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Mock samples were down-sampled and pooled to use as the MACS2 control for each ChIP respectively (see Table S3 in the paper)
   * Peaks that had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
  * chip_noMask_wiDups_tfl1WInputControl_extFragSize_qval10:
   * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads
   * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
   * Inputs from Zhu 2020 were down-sampled and pooled to use as the MACS2 control
   * Peaks that had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
## B. `notebook` directory
jupyter notebooks used to generate some of the plots for the paper:
* ChIPQC_notebook.ipynb - code to parse ChIPQC output
* GSvBL_comparison.ipynb - code to compare lengths of blacklist regions vs green screen regions
* plotting_chip_comparisons.ipynb
 * code to compare the locations of LFY ChIP Seq (Jin 2021) peaks called using various parameter settings (ie. different MACS2 controls, including duplicates (True/False), or artificial masking via data/arabidopsis_blacklist_20inputs.bed blacklist or data/arabidopsis_greenscreen_20inputs.bed green screen) 
  * LFY ChIP-seq (Jin 2021) peaks are compared to:
   * green screen regions (data/arabidopsis_greenscreen_20inputs.bed)
   * LFY ChIP-chip regions (Winter 2011)
  * Possible MACS2 controls:
    * None
    * Input (Jin 2021)
    * Mock (Jin 2021)
 * code to compare the locations of FD ChIP Seq (Zhu 2020) peaks called using various parameter settings (ie. different MACS2 controls, including duplicates (True/False), or artificial masking via data/arabidopsis_blacklist_20inputs.bed blacklist or data/arabidopsis_greenscreen_20inputs.bed green screen) 
  * FD ChIP-seq (Zhu 2020) peaks are compared to:
   * green screen regions (data/arabidopsis_greenscreen_20inputs.bed)
   * FD ChIP-seq regions called using mock control without a mask for artificial regions (Collani 2019)
  * Possible MACS2 controls:
    * None
    * Input (Zhu 2020)
    * Mock (Zhu 2020)
<br>
<br>
# II. Notes to get started
To activate the python environment in this project repository type:
`source env/bin/activate`
