#!/usr/bin/python

import numpy as np
import copy
import heapq
import sys
from scipy.spatial import distance
import csv
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans

if __name__ =='__main__':

    WORK_DIR = '/gpfs/group/su/lhgioia/map/'
    

    points = np.genfromtxt(WORK_DIR + 'results/recount/pca/recount_750_dim_noProj_over50_pc3sd_tpm_log.csv',delimiter=',')
    points = points[:,0:100]

    with open(WORK_DIR + 'results/recount/pca/recount_over50_bulkOnly_pc3sd_filtered_rownames_tpm_log.txt','r') as in_file:
        rownames = in_file.readlines()

    nbrs = NearestNeighbors(n_neighbors=30,algorithm='ball_tree').fit(points)
    distances, indices = nbrs.kneighbors(points)

    with open(WORK_DIR + 'results/recount/clustering/recount_over50_pc3sd_tpm_log_nv100_90th_var_genes_nn_k30.csv','w') as out_file:
        writer = csv.writer(out_file,delimiter=',')
        for i in range(len(points)):
            row = [rownames[i].strip('\n').strip('"')]
            row.extend(indices[i,])
            writer.writerow(row)

    print 'done'
