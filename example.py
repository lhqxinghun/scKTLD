import numpy as np
import sys
import os

import tkinter
import matplotlib
import matplotlib.pyplot as plt
import sklearn.metrics.pairwise as skmetrics
import seaborn as sns

from scKTLD import edge2adj
from scKTLD import callTLD
from scKTLD import displayTLD

## This is for the configuration of color scheme
matplotlib.use('TkAgg')
hiccolors = ["lightyellow", "red"]
my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('hicheat', hiccolors)
matplotlib.cm.register_cmap(cmap = my_cmap)
featurecolors = ["blue", "white", "red"]
my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('featureheat', featurecolors)
matplotlib.cm.register_cmap(cmap = my_cmap)

## Prepare the contact matrix, if the input is sparse, use the following to convert to dense
path_input = "./data/exp-sc/gm12878_cell7_chr3_sparse.txt"
graph_edge = np.loadtxt(path_input)
chr='chr3'
resolution = 50000
graph_adj = edge2adj(graph_edge, chr = chr, resolution = resolution, reference = "hg19")

## call TAD-like domains on the input contact matrix
boundary_spec = callTLD(graph_adj)

## Visualize TAD-like domains on the input or reconstructed contact matrix
displayTLD(graph_adj, boundary_spec, 800, 1000, brecon = True)

plt.show()
#plt.savefig("path to save/filename.tiff", dpi = 350)
#plt.close()
