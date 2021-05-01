This directory contains jupyter notebooks that were used to generate some of the plots for the paper:
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
