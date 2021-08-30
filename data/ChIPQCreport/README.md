This directory contains subdirectories containing results of ChIPQC on inputs and ChIP-seq experiments used in this paper. Below are descriptions of the samples and parameters used to generate the reads and regions that were imported into ChIPQC and output into the specific subdirectory.
* 20inputs_blMask/
  * ChIPQC Results for 20 random inputs (see Table S2 in the paper) after using samtools to remove reads that overlap blacklist (data/arabidopsis_blacklist_20inputs.bed)
* 20inputs_gsMask/
  * ChIPQC Results for 20 random inputs (see Table S2 in the paper) after using samtools to remove reads that overlap greenscreen (data/arabidopsis_greenscreen_20inputs.bed)
* 20inputs_noMask/
  * ChIPQC Results for 20 random inputs (see Table S2 in the paper) 
* chip_blMask_wiDups_bestControlsV5_extFragSize_qval10/ 
  * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads (without duplicates) after blacklist (data/arabidopsis_blacklist_20inputs.bed) masking
  * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
  * Controls from the respective publications were down-sampled and pooled to use as the MACS2 control. If mock and input controls were both provided, input controls were used.
  * Peaks with summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included
* chip_gsMask_wiDups_bestControlsV5_extFragSize_qval10/
  * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads
  * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
  * Controls from the respective publications were down-sampled and pooled to use as the MACS2 control. If mock and input controls were both provided, input controls were used.
  * Peaks that overlap green screen regions (data/arabidopsis_greenscreen_20inputs.bed) OR had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
* chip_noMask_wiDups_bestControlsV5_extFragSize_qval10/
  * This directory contains ChIPQC Results for LFY, TFL1, and FD mapped reads
  * Peaks were called on these reads using MACS2 v2.2.7.1  (—keepdup “auto” —nomodel –extsize [fragment_length] —broad --nolambda  -g 101274395)
  * Controls from the respective publications were down-sampled and pooled to use as the MACS2 control. If mock and input controls were both provided, input controls were used.
  * Peaks that had summit q-values (column 9 in narrowPeak file) greater than 10^-10 were NOT included 
