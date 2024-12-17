
setwd("../data/HER2/BRCA_SpaceRanger")

expr.data <- Seurat::Read10X_h5(filename ="Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5")
obj <- Seurat::CreateSeuratObject(counts = expr.data, project = 'BRCA_HER2', assay = 'Spatial')
obj$slice <- 1
obj$region <- 'BRCA'

#Load the image data

img <- Seurat::Read10X_Image(image.dir = '../data/BRCA_SpaceRanger/spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = obj)]
obj[['image']] <- img

###Standard, dimensionality reduction and clustering

obj <- SCTransform(obj, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)

p <- SpatialDimPlot(obj, label = TRUE, label.size = 3)
saveRDS(obj,file = "../data/HER2/st.obj.rds")