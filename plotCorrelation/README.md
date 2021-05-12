# Dependencies
* pybigwig

# Example commands
To find examples of the format of the input files used for these scripts
see the `example_inputs` directory. Below is how commands using these
example inputs could be formatted:

python3 scripts/coverage_bed_matrix.py \
	[example_inputs/bigwig_paths.csv] \
	[example_inputs/merged_chip_peaks.bed] \
        -o [out_dir] -m [coverage_matrix.csv] 

python3 scripts/readCorrelationPlot.py \
	[out_dir/coverage_matrix.csv] \
	[out_dir/dendrogram_heatmap.png] \
	-lm ward --plot_numbers -k 3 -ri \
	-cf [example_inputs/predictedClusters.csv] \
	-sl [example_inputs/colorLabels.csv]
