
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.spatial import distance_matrix
import novosparc as nc
from scipy import stats
import tangram as tg
from scipy.spatial.distance import cdist
import sys


import os
print(os.getcwd())

import sys
sys.path.append("your_path")
ad_sp = sc.read_h5ad("st.h5ad")
ad_sc = sc.read_h5ad("sc.h5ad")

gene_names = np.array(ad_sc.var.index.values)

dge = ad_sc.to_df().values
num_cells = dge.shape[0]

print ('number of cells and genes in the matrix:', dge.shape)


hvg = np.argsort(np.divide(np.var(dge,axis=0),np.mean(dge,axis=0)+0.0001))
dge_hvg = dge[:,hvg[-2000:]]

try:
    locations = ad_sp.obs[['X','Y']].values
except:
    locations = ad_sp.obs[['x','y']].values
num_locations = locations.shape[0]

p_location, p_expression = nc.rc.create_space_distributions(num_locations, num_cells)
cost_expression, cost_locations = nc.rc.setup_for_OT_reconstruction(dge_hvg,locations,num_neighbors_source = 5,num_neighbors_target = 5)


gene_is=ad_sp.var.index.tolist()
gene_sc=ad_sc.var.index.tolist()
insitu_genes=list(set(gene_is).intersection(gene_sc))

markers_in_sc = np.array([], dtype='int')
for marker in insitu_genes:
    marker_index = np.where(gene_names == marker)[0]
    if len(marker_index) > 0:
        markers_in_sc = np.append(markers_in_sc, marker_index[0])
                
insitu_matrix = np.array(ad_sp.to_df()[insitu_genes])
cost_marker_genes = cdist(dge[:, markers_in_sc]/np.amax(dge[:, markers_in_sc]),insitu_matrix/np.amax(insitu_matrix))

alpha_linear = 0.5
gw = nc.rc._GWadjusted.gromov_wasserstein_adjusted_norm(cost_marker_genes, cost_expression, cost_locations,alpha_linear, p_expression, p_location,'square_loss', epsilon=5e-3, verbose=True)

gamma = gw
for j in range(gamma.shape[1]):
    gamma[:,j] = gamma[:,j]/np.sum(gamma[:,j])
    
ad_map = sc.AnnData(gamma,obs = ad_sc.obs, var=ad_sp.obs)

tg.project_cell_annotations(ad_map, ad_sp, annotation="celltype")

output_file_path="your_path"
ad_sp.obsm['tangram_ct_pred'].to_csv(output_file_path + '/novoSpaRc_decon.csv')


