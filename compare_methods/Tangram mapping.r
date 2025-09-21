###确定性映射
import os

os.chdir("data_path")

import tangram as tg
import scanpy as sc
import numpy as np

ad_sp = sc.read_h5ad("st.h5ad")
ad_sc = sc.read_h5ad("sc.h5ad")

import pandas as pd
rownames = pd.read_csv('rownames.txt', sep='\t')
rownames_first_column = rownames['x'].tolist()

rownames_data = pd.DataFrame({
    'x': rownames_first_column
})
ad_sp.var_names = rownames_data['x']

tg.pp_adatas(ad_sc, ad_sp, genes = None)

ad_map = tg.map_cells_to_space(ad_sc, ad_sp)
ad_map 

ad_map_df = pd.DataFrame(ad_map.X, index=ad_sc.obs.index, columns=ad_sp.obs.index)

ad_map_df.to_csv(output_file_path + '/Tangram_cell.txt')