"""

Cluster the 2D PCA vectors uing kmean clustering

- Please specify the number of clusters


# Example run: cluster PCA
 
singularity exec forgi-centroid-mea-final.simg  python3  draw-kmeans-clusters.py --input_file COOLAIR3_R1_R2-forgi-pca.csv --tag  COOLAIR3_R1_R2 --num_clusters 3

"""

import matplotlib as mpl
mpl.use('Agg')

import sys
import os
import argparse
import jitu
from collections import defaultdict
import subprocess

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min


def handler():
    """
    Get command line inputs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file",   type=str,  required=True)
    parser.add_argument("-t", "--tag",          type=str, required=True)
    parser.add_argument("-n", "--num_clusters", type=int, default=3)


    return parser.parse_args()


if __name__ == "__main__":

    args = handler()

    # number of clusters
    num_clusters = args.num_clusters or 3

    # load data from file     
    df = pd.read_csv(args.input_file, index_col = "structure")
    print ( df.head())

    # on 2D-PCA do kmean to colour points
    kmeans = KMeans(n_clusters = num_clusters, random_state=0).fit(df)
 
    centroids =  kmeans.cluster_centers_.tolist()
    k_clust_ids = kmeans.labels_.tolist()

    print( 'using PCA coordinates:')
    for i in range(num_clusters):
       print(f'number_of_points_in_cluster{i+1} = ', k_clust_ids.count(i) )
    print('')
    
    closest, _ = pairwise_distances_argmin_min( kmeans.cluster_centers_, df)

    # colour for each cluster
    cluster_color = { 0: '#4daf4a', 
                      1: '#ff7f00', 
                      2: '#f781bf'} # greenish, orange, pinkish

    # check human colour https://www.w3schools.com/colors/colors_hexadecimal.asp
    human_color =  {'#4daf4a': 'green', 
                    '#ff7f00': 'orange', 
                    '#f781bf': 'pink' 
    }

    print( 'Most representative member using PCA coordinates:')
    for k, i in enumerate(closest):
        print ('cluster' + str(k+1), 'position_zero_based:' + str(i),  'structure_id:', df.index[i], 'colour: ', cluster_color[k], human_color[cluster_color[k]] )
     
    # plotting
    fig = plt.figure(figsize = (6,6))
    ax = fig.add_subplot(1,1,1) 

    plt.scatter(df['PC1'], df['PC2'], c = [cluster_color[x] for x in kmeans.labels_], s = 50, alpha = 0.50 )#, alpha=0.5) # 1: yellow 0;purple   colormap='viridis'

    for k, i in enumerate(closest): 
        # add red circle    
        plt.scatter(df.iloc[i]['PC1'], df.iloc[i]['PC2'],  c =  cluster_color[k], s = 50, edgecolors= '#e41a1c')#, mfc = 'r')

        # add text annotation
        plt.text(df.iloc[i]['PC1'] + 1.5, df.iloc[i]['PC2']  + 1.5, df.index[i], horizontalalignment='left',  color='black',   fontsize=8 )

    ax.set_xlabel('PC 1', fontsize = 15)
    ax.set_ylabel('PC 2', fontsize = 15)
    ax.set_title('Kmean clustering on PCA forgi vectors', fontsize = 15)
    ax.set_aspect('equal', 'box')
    plt.axis('equal')

    plt.savefig(args.tag + '-clusters.png', dpi=1200)
    plt.savefig(args.tag + '-clusters.pdf')

    # output to a csv file
    df['cluster_id'] =  kmeans.labels_.tolist()
    df['cluster_name'] = [ str(k+1)   for k in kmeans.labels_.tolist()]
    df['cluster_colour'] = [ human_color[cluster_color[k]]  for k in kmeans.labels_.tolist()]
    df['cluster_hex_colour'] = [ cluster_color[k]  for k in kmeans.labels_.tolist()]

    df.to_csv(args.tag + '-clusters.csv', index = True)

    print ('Done')

"""

[cheemaj@NBI-HPC interactive test-sm-structure-seq-pipeline-run]$ singularity exec forgi-centroid-mea-final.simg  python3  mark-pca-groups.py --input_file COOLAIR3_R1_R2-forgi-pca.csv --tag  COOLAIR3_R1_R2
                 PC1       PC2
structure
bit_1       9.224530  8.650775
bit_2      20.918670  2.865675
bit_3       6.747810  5.819560
bit_4       9.161598  8.607520
bit_5      16.818061 -0.676877

using PCA coordinates:
 number_of_points_in_cluster1 =  152
 number_of_points_in_cluster2 =  460
 number_of_points_in_cluster3 =  190

Most representative using PCA coordinates:

cluster1 position_zero_based:8 structure_id: PC1     0.387908
 PC2    32.251042
 Name: bit_9, dtype: float64 bit_9

cluster2 position_zero_based:263 structure_id: PC1   -23.765013
 PC2    -8.086184
 Name: bit_264, dtype: float64 bit_264
 cluster3 position_zero_based:157 structure_id: PC1    53.528874
 PC2    -5.985844
 Name: bit_158, dtype: float64 bit_158

Done
[cheemaj@NBI-HPC interactive test-sm-structure-seq-pipeline-run]$

"""

