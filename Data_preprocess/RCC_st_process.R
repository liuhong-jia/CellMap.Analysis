####RCC spatial transcriptomics data processing

expr.data <- Seurat::Read10X_h5(filename ="../data/RCC/GSM5924052_frozen_c_3_filtered_feature_bc_matrix.h5")
obj <- Seurat::CreateSeuratObject(counts = expr.data, project = 'RCC', assay = 'Spatial')
obj$slice <- 1
obj$region <- 'BRCA'

#Load the image data

img <- Seurat::Read10X_Image(image.dir = 'spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = obj)]
obj[['image']] <- img

TLS_data <- read.csv("../data/TLS/TLS_annotation.csv")
head(st.obj@meta.data)
st.obj$TLS <- TLS_data$TLS_2_cat

###Standard, dimensionality reduction and clustering

obj <- SCTransform(obj, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)



saveRDS(obj,file = "../data/RCC/st.obj.rds")