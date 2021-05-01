This directory contains subdirectories containing results of ChIPQC on inputs and ChIP-seq experiments used in this paper. Below are descriptions of the samples and parameters used to generate the reads and regions that were imported into ChIPQC and output into the specific subdirectory.
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
