

library(Seurat)
library(SeuratDisk)
sc <- readRDS("your_path/sc.obj.rds")
st <- readRDS("your_path/st.obj.rds")

options(Seurat.object.assay.version = "v5")
sc = CreateSeuratObject(counts = GetAssayData(sc, slot = "counts"),
		                                meta.data = sc@meta.data)

st = CreateSeuratObject(counts = GetAssayData(st, slot = "counts"),
		                                meta.data = st@meta.data)


SaveH5Seurat(sc,filename = "sc.h5Seurat")
Convert("sc.h5Seurat", dest = "h5ad")

SaveH5Seurat(st,filename = "st.h5Seurat")
Convert("st.h5Seurat", dest = "h5ad")