#!/usr/bin/env python3

# -*- coding: iso-8859-15 -*-
# 2019, Samantha Klasfeld, the Wagner Lab
# the Perelman School of Medicine, the University of Pennsylvania
# Samantha Klasfeld, 12-9-2019

import os
import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.patches as mpatches



parser = argparse.ArgumentParser(description="Using a counts matrix \
    for each sample in specific regions compare each sample \
    using agglomerative (bottom-up) hierarchical clustering")
parser.add_argument('infile', help='a comma-delimited \
    counts matrix formatted with columns `chrom`,`start`,`stop` \
    and sample names for each sample \
    (eg.chrom,start,stop,samp1,samp2,samp3).')
parser.add_argument('out', help='file name of output plot')
parser.add_argument('-cm','--corr_method', help='method used \
    to calculate pairwise correlation (default= "pearson")', choices=['pearson', \
    'kendall', 'spearman'], default="pearson")
parser.add_argument('-lm','--link_method', help='linkage \
    methods for calculating the distance between the \
    newly formed cluster and all the other clusters. See \
    scipy.cluster.hierarchy.linkage documentation for more \
    information', choices=['single', 'complete', 'average', \
    'weighted', 'centroid', 'median', 'ward'], default="complete")
parser.add_argument('--plot_numbers', help='plot numbers on \
    heatmap', action='store_true', default=False)
parser.add_argument('-k','--k_clusters', help='number of clusters to color \
    (default:2)', type=int, default=2)
parser.add_argument('-cf','--cluster_file',type=str,
    help='set a file-path that contains a table with cluster \
    information of the samples.The file should include \
    comma-delimited columns. The first column labeled "sample_name" \
    should include all of the sample names that match the infiles. \
    The rest of the columns should include labels starting from 1 \
    for the clusters you want to compare. Name each clustering \
    method in the header.')
parser.add_argument('-ri','--calcRandIndexForSpecificClustering', 
    action="store_true",
    help = 'Set this if you want to calculate the rand-index for \
    a given set of clusters. If this is NOT set, we will print the \
    RandIndex comparing the two clustering algorithms.')
parser.add_argument('-e','--extraInfo', type=str,
    help = 'If `calcRandIndexForSpecificClustering` is set, \
    the user may want to run a batch script to make multiple \
    measurements. Set this to a label to identify the params \
    used in this code. For example, you may want to print out \
    the params used to make a blacklist/greenscreen that was \
    used in the samples.')
parser.add_argument('-sl','--samplabels',
    help = 'set a file-path that contains a table with sample \
    meta data for the samples.The file should include \
    comma-delimited columns. The first column labeled "sample_name" \
    should include all of the sample names that match the infiles. \
    The rest of the columns should include discrete meta-data fields \
    to identify the sample. Name columns `color` and/or `shape` \
    respectively to label a specific sample by that color and/or \
    marker shape. See https://matplotlib.org/3.1.1/api/markers_api.html \
    for marker syntax. Warning: cannot use the pixel marker.')
parser.add_argument('--vmin_max',nargs=2,default=[None,None],
    help = 'The colorbar range for heatmap')
args = parser.parse_args()

colormap='rainbow'

if args.calcRandIndexForSpecificClustering:
    if os.path.isfile(args.cluster_file):
        trueClusters_df = pd.read_csv(
            args.cluster_file, 
            sep=",", index_col=0)
    else:
        sys.exit("ERROR: cannot find file: %s" % 
            args.cluster_file)

if args.samplabels:
    if os.path.isfile(args.samplabels):
        samplabels_df = pd.read_csv(
            args.samplabels, 
            sep=",")
    else:
        sys.exit("ERROR: cannot find file: %s" % 
            args.samplabels)

# create and export figure
fig = plt.figure(figsize=(11, 9.5))

counts_df = pd.read_csv(args.infile, sep=",")
region_cols = ["chrom","start","stop"]
for rc in region_cols:
    if rc not in list(counts_df.columns):
        sys.exit(("ERROR: Matrix in %s does"+
            " not contain the %s column.") % 
        (args.infile, rc))
counts_df = counts_df.astype(
    {'start': 'int64','stop': 'int64'})
counts_df = counts_df.set_index(region_cols)
corr_matrix = counts_df.corr(method =args.corr_method)
edge_color = 'black'
if corr_matrix.shape[0] > 30:
    # when there are too many rows it is better to remove
    # the black lines surrounding the boxes in the heatmap
    edge_color = 'none'

# get color map
cmap = plt.get_cmap(colormap)

# initialize axes
if not args.samplabels:
    dend_ax = plt.subplot2grid((10, 10), (0, 0), colspan=1, rowspan=9)
    heatmap_ax = plt.subplot2grid((10, 10), (0, 1), colspan=9, rowspan=9)
    colbar_ax = plt.subplot2grid((10, 10), (9, 0), colspan=10, rowspan=1)
else:
    label_ax1 = plt.subplot2grid((10, 10), (0, 1), colspan=8, rowspan=1)
    dend_ax = plt.subplot2grid((10, 10), (1, 0), colspan=1, rowspan=8)
    heatmap_ax = plt.subplot2grid((10, 10), (1, 1), colspan=8, rowspan=8)
    label_ax2 = plt.subplot2grid((10,10), (1,9), colspan=1, rowspan=8)
    colbar_ax = plt.subplot2grid((10, 10), (9, 0), colspan=9, rowspan=1)    

# calculate linkage
linked = linkage(corr_matrix, args.link_method)

color_fclust = fcluster(linked,args.k_clusters,"maxclust")

# calculate labels for dendrogram
sampList = list(corr_matrix.columns)
n=len(sampList)
labels=list('' for i in range(n))
for i in range(n):
    labels[i]=str(i)+ ',' + str(color_fclust[i])

# calculate color threshold for dendrogram
ct=linked[-(args.k_clusters-1),2]


for k in range(2,n+1):
    fclust = fcluster(linked,k,"maxclust")
    colname="0_"+str(k)
    if k==2:
        cluster_df = pd.DataFrame({"sample_name":sampList,
            colname:fclust})
    else:
        cluster_df[colname]=fclust
cluster_df.set_index('sample_name', inplace=True)




# plot dendrogram
d_plot = dendrogram(linked, 
    orientation='left',
    distance_sort='ascending', 
    ax=dend_ax,
    labels=labels,color_threshold=ct)
# remove grid behind dendrogram
dend_ax.set_axis_off()


if args.plot_numbers:
    cmap = cmap.from_list(colormap + "clipped",
        cmap(np.linspace(0, 0.9, 10)))

# reorder matrix to match dendrogram
index = d_plot['leaves']


corr_matrix_denOrdered = corr_matrix.iloc[index, :]
corr_matrix_denOrdered = corr_matrix_denOrdered.iloc[:, index]

if args.samplabels:
    denOrder_dic={}
    for count, value in enumerate(list(corr_matrix_denOrdered.columns)):
        denOrder_dic[value]=count
    samplabels_df["s_idx"] = samplabels_df["sample_name"].apply(lambda x:denOrder_dic[x])
    samplabels_df = samplabels_df.set_index("s_idx")
    samplabels_df = samplabels_df.sort_index()
    samplabels_df = samplabels_df.set_index("sample_name")
# plot heatmap
img_mat = heatmap_ax.pcolormesh(corr_matrix_denOrdered,
    cmap=cmap,vmin=args.vmin_max[0],
    vmax=args.vmin_max[1])
## turn off default y-ticks
heatmap_ax.set_yticks([])

# set labels to heatmap
if not args.samplabels:
    heatmap_ax.yaxis.tick_right()
    heatmap_ax.set_yticks(np.arange(corr_matrix .shape[0]) + 0.5)
    heatmap_ax.set_yticklabels(np.array(counts_df.columns)[index], 
                fontsize='x-large')

    heatmap_ax.xaxis.set_tick_params(labeltop=True)
    heatmap_ax.xaxis.set_tick_params(labelbottom=False)
    heatmap_ax.set_xticks(np.arange(corr_matrix .shape[0]) + 0.5)
    heatmap_ax.set_xticklabels(np.array(counts_df.columns)[index], 
        rotation=45, ha='left')
else:
    sns.set(font_scale=3)
    if not ("color" in list(samplabels_df.columns) or 
        "shape" in list(samplabels_df.columns)):
        # create dictionary with value to integer mappings
        value_to_int = {value: i for i, value in enumerate(sorted(pd.unique(samplabels_df.values.ravel())))}
        hm = sns.heatmap(samplabels_df.replace(value_to_int).T, 
            cmap=plt.get_cmap('jet'), ax=label_ax1, cbar=False, linewidths=.5)
        rev_samplabels_df = samplabels_df[samplabels_df.columns[::-1]]
        hm2 = sns.heatmap(rev_samplabels_df.replace(value_to_int).reset_index(drop=True).sort_index(ascending=False, axis=0), 
            cmap=plt.get_cmap('jet'), ax=label_ax2, cbar=False, linewidths=.5)
        # hide x labels
        hm.set(xticklabels=[])
        hm.tick_params(bottom=False,labelsize=25)
        hm2.set(yticklabels=[])
        hm2.tick_params(left=False,labelsize=15,labelrotation=70)
        # add legend
        box = label_ax2.get_position()
        label_ax2.set_position([box.x0, box.y0, box.width * 0.7, box.height])
        legend_ax = fig.add_axes([1, .5, 1, .1])
        legend_ax.axis('off')
        # reconstruct color map
        colors = plt.cm.jet(np.linspace(0, 1, len(value_to_int)))
        # add color map to legend
        patches = [mpatches.Patch(facecolor=c, edgecolor=c) for c in colors]
        legend = legend_ax.legend(patches,
            sorted(value_to_int.keys()),
            handlelength=0.8, loc='lower left')
        for t in legend.get_texts():
            t.set_ha("left")
    else:
        if "color" in list(samplabels_df.columns):
            color_list=list(samplabels_df["color"])
            samplabels_df = samplabels_df.drop(["color"], axis=1)
        else:
            color_list=["black" for x in range(0,len(samplabels_df))]
        if  "shape" in list(samplabels_df.columns):
            shape_list=list(samplabels_df["shape"])
            samplabels_df = samplabels_df.drop(["shape"], axis=1)
        else:
            shape_list=["s" for x in range(0,len(samplabels_df))]
        
        label_list = list(samplabels_df.apply(lambda x: ",".join(x), axis=1))

        for i in range(0,len(samplabels_df)):
            y_tick = [1]
            x_tick = [i+.5]
        
            label_ax1.plot(x_tick, y_tick,
                marker=shape_list[i], 
                c=color_list[i],
                label=label_list[i],
                markersize=12)
            label_ax1.set_xlim(left=0,right=len(samplabels_df))
            label_ax1.set_facecolor('white')
            label_ax1.axis('off')
            label_ax2.plot(y_tick, x_tick,
                marker=shape_list[i], 
                c=color_list[i],
                label=label_list[i],
                markersize=12)
            label_ax2.set_ylim(bottom=0,top=len(samplabels_df))
            label_ax2.set_facecolor('white')
            label_ax2.axis('off')

# Plot colorbar
cobar = plt.colorbar(img_mat, cax=colbar_ax, orientation='horizontal')
cobar.solids.set_edgecolor("face")

# Plot Numbers
if args.plot_numbers:
    num_rows=len(list(counts_df.columns))
    # get correct font size
    # set a font size according to figure length
    if num_rows < 6:
        font_size = 14
    elif num_rows > 40:
        font_size = 5
    else:
        font_size = int(14 - 0.25 * num_rows)
   # matplotlib.rcParams.update({'font.size': font_size})
    # write numbers    
    for row in range(num_rows):
        for col in range(num_rows):
            heatmap_ax.text(row + 0.5, col + 0.5,
                          "{:.2f}".format(corr_matrix_denOrdered.iloc[row, col]),
                          ha='center', va='center',
                          fontsize=font_size)
    heatmap_ax.axes.get_xaxis().set_visible(False)

plt.savefig(args.out, bbox_inches='tight')

if args.calcRandIndexForSpecificClustering:
    if args.extraInfo:
        sys.stdout.write("> meta: %s\n" % (args.extraInfo))
    cluster_types = list(trueClusters_df.columns)
    for clustType in cluster_types:
        k = np.max(trueClusters_df[clustType])
        # match_count = count of pairs of elements that are in the 
        # same subset in clustMethod1 and in the same subset in clustMethod2
        # plus the number of pairs of elements that are in 
        # different subsets in clustMethod1 and in different subsets in 
        # clustMethod2
        match_count = 0.0
        total_count = 0.0
        for samp1 in sampList:
            for samp2 in sampList:
                # check if both are in the same subset
                sameClusterInGeneratedClusters = \
                    (cluster_df.loc[samp1, "0_"+str(k)] == 
                    cluster_df.loc[samp2, "0_"+str(k)])
                sameClusterInTrueClusters = \
                    (trueClusters_df.loc[samp1, clustType] == 
                        trueClusters_df.loc[samp2, clustType])

                if (sameClusterInGeneratedClusters == 
                    sameClusterInTrueClusters):
                    match_count = match_count + 1.0
                total_count = total_count + 1.0
        randIndex_value = ( (match_count) / (total_count) )
        sys.stdout.write("> %s: %.2f\n" % (clustType, randIndex_value))
        match_count=0.0
        total_count=0.0
