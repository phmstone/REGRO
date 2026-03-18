####################################################################################################
# This script:
# 1. Imports presence/absence .TSV file made by presenceAbsence.py
# 2. Uses the numeric values to make a "heatmap" figure
# 3. Saves the heatmap figure according to name given by user on command line or default name
####################################################################################################

# ----------------------------------------------------------
# Import packages
# ----------------------------------------------------------

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.colorbar import ColorbarBase
import matplotlib.patches as mpatches


# ----------------------------------------------------------
# Command line arguments
# ----------------------------------------------------------

parser = argparse.ArgumentParser(description="Plot gene status matrix from TSV file.")
parser.add_argument("--input", help="Input presence/absence TSV file ")
parser.add_argument("--output", default="heatMapPlot.png", help="Output image filename (default: heatMapPlot.png)")
# parser.add_argument("--removeRows", help="GenBank ID for accessions in TSV that should be left out of the plot")

args = parser.parse_args()

# ----------------------------------------------------------
# Load presence/absence .TSV
# ----------------------------------------------------------

# generated with presenceAbsence.py

# read in .TSV file as a dataframe
df = pd.read_csv(args.input, sep="\t")

# combine taxon and GenBank ID into a new row name for the plot
taxon_col = df.columns[0]
genbank_col = df.columns[1]
row_labels = df[taxon_col] + " (" + df[genbank_col] + ")"

# extract gene presence/absence info data (starts in column 3)
gene_data = df.iloc[:, 2:]
gene_data.index = row_labels


# ----------------------------------------------------------
# Sort out colours
# ----------------------------------------------------------

# at the moment 0 is present and 1 is absent, might switch these around
code_labels = {
    0: "Present",
    1: "Absent",
    2: "Pseudogene"}

# colour scheme 
colors = [
    "white",   # 0 Present
    "black",   # 1 Absent
    "gray"]    # 2 Pseudogene

# make the colour scheme into a colour map for matplotlib
cmap = ListedColormap(colors)
norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

# ----------------------------------------------------------
# Make the heatmap
# ----------------------------------------------------------

# make the figure the size of number of rows and columns x 0.6 inches
plt.figure(figsize=(len(gene_data.columns) * 0.6, 
                    len(gene_data.index) * 0.6))

# make a heatmap with seaborn
heatmap = sns.heatmap(
    gene_data, 			# the 1s 0s and 2s
    cmap=cmap,			# using assigned colour map
    norm=norm,			# map to specific colours
    linewidths=0.5,		# thin grid lines
    linecolor="black",	# black grid lines
    cbar=False,			# no color bar!!!
    square=True) 		# square grid cells

# ----------------------------------------------------------
# Make the legend
# ----------------------------------------------------------

# make discrete legend patches for each colour
patches = [
    mpatches.Patch(color=colors[0], label=code_labels[0]),
    mpatches.Patch(color=colors[1], label=code_labels[1]),
    mpatches.Patch(color=colors[2], label=code_labels[2])]

# add the legend outside the heatmap
ax = heatmap  #  matplotlib Axes
ax.legend(
    handles=patches,
    title="Gene Status",
    loc="center left",
    bbox_to_anchor=(1, 0.5), # position the legend halfway down the plot and to the right
    frameon=True,
    fontsize=18,
    title_fontsize=24)


# ----------------------------------------------------------
# Labels and plot format
# ----------------------------------------------------------

# plot labels
plt.xlabel("Genes", fontsize=24)
plt.ylabel("Taxon (GenBank ID)", fontsize=24)
plt.title("Gene Status by Taxon", fontsize=36)

# axes labels
heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=18) # rotate x axis labels to be able to read gene names
heatmap.set_yticklabels(heatmap.get_yticklabels(), fontsize=18)

plt.tight_layout() # minimise the amount of white space in the plot

# ----------------------------------------------------------
# Save the plot
# ----------------------------------------------------------

plt.savefig(args.output, dpi=300, bbox_inches="tight") # save plot according to given name, at 300 dpi resolution, and reducing whitespace
plt.close()

print(f"Plot saved as {args.output}")
