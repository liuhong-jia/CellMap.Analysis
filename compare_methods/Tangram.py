import tangram as tg
import scanpy as sc
import numpy as np

ad_sp = sc.read_h5ad("st.h5ad")
ad_sc = sc.read_h5ad("sc.h5ad")

tg.pp_adatas(ad_sc, ad_sp, genes = None)
ad_map = tg.map_cells_to_space(ad_sc, ad_sp)
spotID = ad_map.var.index[np.argmax(ad_map.X, axis = 1)]
ad_map.obs['SpotID'] = spotID

meta_data = ad_map.obs
meta_data['id'] = meta_data.index

subset_data = meta_data[["celltype","id", "SpotID"]]
subset_data.to_csv("your_path/tangram.csv",index = True,sep = ",")