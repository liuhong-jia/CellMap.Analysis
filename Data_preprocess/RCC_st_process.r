####RCC spatial transcriptomics data processing

expr.data <- Seurat::Read10X_h5(filename ="GSM5924052_frozen_c_3_filtered_feature_bc_matrix.h5")
obj <- Seurat::CreateSeuratObject(counts = expr.data, project = 'RCC', assay = 'Spatial')
obj$slice <- 1
obj$region <- 'BRCA'

#Load the image data

img <- Seurat::Read10X_Image(image.dir = 'spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = obj)]
obj[['image']] <- img

###Standard, dimensionality reduction and clustering

obj <- SCTransform(obj, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)

saveRDS(obj,file = "st.obj.rds")