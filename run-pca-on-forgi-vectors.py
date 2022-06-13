"""
Run PCA 

singularity exec forgi-centroid-mea-final.simg  python3  run-pca-on-forgi-vectors.py --input_file forgi-vect-ser.txt --tag COOLAIR3_R1_R2 --csv_file COOLAIR3_R1_R2-forgi-pca.csv

"""

import sys
import matplotlib as mpl
mpl.use('Agg')

import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from collections import defaultdict


def handler():
    """
    Get command line inputs
    """ 
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=argparse.FileType('r'),  required=True)
    parser.add_argument("-t", "--tag",  type=str, required=True)
    parser.add_argument("-c", "--csv_file",  type=str, required=True)

    # 'forgi-vect-ser.txt'
    #parser.set_defaults(input_file = 'forgi-vect-ser.txt',
    #                    csv_file = 'COOLAIR3_R1_R2-forgi-pca.csv',
    #                    tag = 'COOLAIR3_R1_R2')

    args = parser.parse_args()

    return args


if __name__ == "__main__":


   args = handler()

   #-----------------
   data =[]
   ser2head = {}
   headers = []
   sizerForgi = defaultdict(int)
   
   for i, line in enumerate(args.input_file):
       head, digits = line.strip().split() 
       A = [float(d) for d in  list(digits)]
       data.append(A)
       ser2head[i] = head
       headers.append(head) 
       sizerForgi[digits] += 1

   print ('Input_forgi_vectors:', len(ser2head))
   print ('uniqForgivectors:', len(sizerForgi))

   #--------------- 
   x = np.asarray(data, dtype='float64')

   pca = PCA(n_components=2)#
   # StandardScaler()
   principalComponents = pca.fit_transform(x)

   print ('Explained_variance_ratio_:', pca.explained_variance_ratio_.tolist())
   print ('Explained_variance_ratio_cumsum:', pca.explained_variance_ratio_.cumsum().tolist())

   # write PCA to a csv file 
   principalDf = pd.DataFrame(data = principalComponents, columns = ['PC1', 'PC2'])
   principalDf['structure'] =  headers
   principalDf.to_csv(args.csv_file, index = False)
 

   # plotting
   fig = plt.figure(figsize = (8,8))
   ax = fig.add_subplot(1,1,1) 

   plt.scatter(principalComponents[:, 0], principalComponents[:, 1], c='#ff7f00', s = 10.0 )


   ax.set_xlabel('Principal Component 1', fontsize = 15)
   ax.set_ylabel('Principal Component 2', fontsize = 15)
   ax.set_title('PCA on contrafolded default Forgi vectors', fontsize = 15)


   plt.savefig(args.tag + '.png', dpi = 600)
   plt.savefig(args.tag + '.pdf')

   print ('DONE')
