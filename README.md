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
  * Twenty input-seq samples from *Arabidopsis thaliana* were cleaned and mapped (MAPQ>=30). See Table S2 in the Greenscreen Paper.
* data/arabidopsis_greenscreen_20inputs.bed
 * Green screen regions were generated using the green screen pipeline.
   1. MACS2 v2.2.7.1 (—keepdup “auto” —no model –extsize [read length] —broad --nolambda  -g 101274395) was run on each of the twenty inputs (Table S2 in the paper). 
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
* ChIPQCreport/ : directory containing the results of ChIPQC on inputs and ChIP-seq experiments used in this paper
## B. `notebook` directory
jupyter notebooks used to generate some of the plots for the paper 
<br>
<br>
# II. Notes to get started
To activate the python environment in this project repository type:
`source env/bin/activate`
