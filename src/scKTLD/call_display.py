import numpy as np
import sys
import os

import matplotlib
import matplotlib.pyplot as plt
import sklearn.metrics.pairwise as skmetrics
import seaborn as sns

#sys.path.append("./")   #To debug this module, change the work directory to the path which contains scKTLD
from scKTLD.preprocess import edge2adj
from scKTLD.nodeEmb import ProNE
from scKTLD.cpDector import cpDector

step = 10
theta = 0.5
mu = 0.2

def callTLD(graph, dimension = 128, penalty = 1.8, brecon=False):

    model = ProNE(graph,  dimension)
    features_matrix = model.pre_factorization(model.matrix0, model.matrix0)
    embeddings_matrix = model.chebyshev_gaussian(model.matrix0, features_matrix, step, mu, theta)
    boundary_spec = cpDector(embeddings_matrix, pen = penalty)
    if brecon:
        recon = skmetrics.linear_kernel(embeddings_matrix)
        recon = (recon-np.min(recon))/(np.max(recon)-np.min(recon))
        return boundary_spec, recon
    else:
        return boundary_spec

def displayTLD(graph, boundary_spec, start, stop, brecon=True, dimension = 128):

    ### define the color scheme
    hiccolors = ["lightyellow", "red"]
    my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('hicheat', hiccolors)
    matplotlib.cm.register_cmap(cmap = my_cmap)
    featurecolors = ["blue", "white", "red"]
    my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('featureheat', featurecolors)
    matplotlib.cm.register_cmap(cmap = my_cmap)

    ### graph segmentation
    graph_seg = graph[start:stop, start:stop]

    ### boundary segmentation
    if (boundary_spec is not None) & (len(boundary_spec)!=0):
        bd_temp_ds= boundary_spec - start
        bd_temp_ds = bd_temp_ds[np.bitwise_and(bd_temp_ds >= 0, bd_temp_ds < (stop-start))]
        bd_temp_ds = np.append(np.insert(bd_temp_ds,0,0),stop-start-1)
        boundary_plot_ds = np.unique(bd_temp_ds)+0.5
    else:
        boundary_plot_ds = []

    ### visualization
    if brecon:
        model = ProNE(graph,  dimension)
        features_matrix = model.pre_factorization(model.matrix0, model.matrix0)
        embeddings_matrix = model.chebyshev_gaussian(model.matrix0, features_matrix, step, mu, theta)
        recon = skmetrics.linear_kernel(embeddings_matrix)
        recon = (recon-np.min(recon))/(np.max(recon)-np.min(recon))
        recon_seg = recon[start:stop, start:stop]

        fig = plt.figure(figsize=(13.0/2.54,5.0/2.54), constrained_layout=False, dpi=150)
        plt.rcParams['font.sans-serif']='Liberation Sans'
        plt.rcParams['font.size']=8
        #
        plt.rcParams['axes.xmargin']=0.01
        plt.rcParams['axes.xmargin']=0.01
        plt.rcParams['savefig.pad_inches']=0.01
        plt.rcParams['savefig.bbox']='tight'

        plt.subplot2grid((1,2),(0,0), rowspan = 1, colspan = 1)
        sns.heatmap(np.log10(graph_seg+1), cmap='hicheat', xticklabels = False, yticklabels = False)
        currentAxis = plt.gca()
        for i in range(len(boundary_plot_ds)-1):
            rect = plt.Rectangle((boundary_plot_ds[i], boundary_plot_ds[i]), boundary_plot_ds[i+1]-boundary_plot_ds[i],boundary_plot_ds[i+1]-boundary_plot_ds[i], linewidth = 1, edgecolor = 'b', facecolor = "none")
            currentAxis.add_patch(rect)

        plt.subplot2grid((1,2),(0,1), rowspan = 1, colspan = 1)
        sns.heatmap(recon_seg, cmap='hicheat', xticklabels = False, yticklabels = False)
        currentAxis = plt.gca()
        for i in range(len(boundary_plot_ds)-1):
            rect = plt.Rectangle((boundary_plot_ds[i], boundary_plot_ds[i]), boundary_plot_ds[i+1]-boundary_plot_ds[i],boundary_plot_ds[i+1]-boundary_plot_ds[i], linewidth = 1, edgecolor = 'b', facecolor = "none")
            currentAxis.add_patch(rect)
    else:
        fig = plt.figure(figsize=(6.5/2.54,5.0/2.54), constrained_layout=False, dpi=150)
        plt.rcParams['font.sans-serif']='Liberation Sans'
        plt.rcParams['font.size']=8
        #
        plt.rcParams['axes.xmargin']=0.01
        plt.rcParams['axes.xmargin']=0.01
        plt.rcParams['savefig.pad_inches']=0.01
        plt.rcParams['savefig.bbox']='tight'
        sns.heatmap(np.log10(graph_seg+1), cmap='hicheat', xticklabels = False, yticklabels = False)
        currentAxis = plt.gca()
        for i in range(len(boundary_plot_ds)-1):
            rect = plt.Rectangle((boundary_plot_ds[i], boundary_plot_ds[i]), boundary_plot_ds[i+1]-boundary_plot_ds[i],boundary_plot_ds[i+1]-boundary_plot_ds[i], linewidth = 1, edgecolor = 'b', facecolor = "none")
            currentAxis.add_patch(rect)

if __name__ == '__main__':
     print(callTLD)
