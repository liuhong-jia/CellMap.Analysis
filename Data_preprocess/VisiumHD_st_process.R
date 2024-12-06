###Visium HD st data preprocess
library(Seurat)
library(hdf5r)
st.data <- Load10X_Spatial(data.dir = getwd(), bin.size = 8)
DefaultAssay(st.data) <- "Spatial.008um"

st.data <- NormalizeData(st.data)
st.data <- FindVariableFeatures(st.data)
st.data <- ScaleData(st.data)

#object <- SCTransform(object, assay = "Spatial.008um", return.only.var.genes = FALSE, verbose = FALSE)
st.data <- RunPCA(st.data, assay = "Spatial.008um", reduction.name = "pca")
st.data <- FindNeighbors(st.data, assay = "Spatial.008um", reduction = "pca", dims = 1:30)
st.data <- FindClusters(st.data, cluster.name = "seurat_cluster.008um", resolution = 0.8)



metadata <- read.csv("DeconvolutionResults_P1CRC.csv")

rownames(metadata) <- metadata$barcode
metadata <- metadata[rownames(st.data@meta.data),]
st.data <- AddMetaData(st.data, metadata = metadata)
st.data <- subset(st.data, subset = DeconvolutionClass == "singlet")

counts <- GetAssayData(st.data,layer = "counts")
coord.df <- st.data@images$slice1.008um$centroids@coords %>% as.data.frame
rownames(coord.df) <- st.data@images$slice1.008um@boundaries$centroids@cells
metadata = st.data@meta.data

createSpObj <- function(counts, coord.df, coord.label = c("x", "y"), meta.data = NULL,class = "SlideSeq") {
    matched_spots <- intersect(colnames(counts), rownames(coord.df))
    if (!is.null(meta.data)) {
        matched_spots <- intersect(matched_spots, rownames(meta.data))
    }
    if (length(matched_spots) == 0) {
        stop("No matching barcodes found between the count matrix and the coordinates.")
    }

    obj <- CreateSeuratObject(counts = counts[, matched_spots, drop = FALSE], assay = "Spatial")
    obj@images$image <- new(
	    Class = class,
        assay = "Spatial",
        key = "image_",
        coordinates = coord.df[matched_spots, coord.label, drop = FALSE]
    )
    if (!is.null(meta.data)) {
        obj@meta.data <- cbind(obj@meta.data, meta.data[matched_spots, , drop = FALSE])
    }
    return(obj)
}

st.obj <- createSpObj(counts, coord.df, coord.label = c("x", "y"), meta.data = metadata)

library(tidydr)
library(Seurat)
library(dplyr)

downSamplSeurat <- function(obj, cluster.col = NULL, cnt = 200, seed = 123, percent = 0.3) {
  set.seed(seed)
  if (!is.null(cluster.col)) Idents(obj) <- obj@meta.data[, cluster.col]
  cells <- Idents(obj) %>% table
  sub.cells <- sapply(names(cells), function(xx) {
    sub.cells <- Idents(obj)[Idents(obj) == xx] %>% names
    cnt <- ifelse(is.null(percent), cnt, length(sub.cells) * percent)
    if (length(sub.cells) > cnt) sub.cells <- sample(sub.cells, cnt, replace = FALSE)
    return(sub.cells)
  }) %>% unlist(use.names = F)
  subset(obj, cells = sub.cells)
}

library(dplyr)

obj2 <- downSamplSeurat(st.obj,cluster.col = NULL,cnt = 2000,seed = 123,percent = NULL)
obj2 <- NormalizeData(obj2, assay = "Spatial", verbose = FALSE)
obj2 <- FindVariableFeatures(obj2, assay = "Spatial", selection.method = "vst", nfeatures = 2000)
obj2 <- ScaleData(obj2,assay = "Spatial") 
obj2 <- RunPCA(obj2,verbose = FALSE)
obj2 <- FindNeighbors(obj2, reduction = "pca", dims = 1:30)
obj2 <- FindClusters(obj2, verbose = FALSE,resolution = 0.5)
obj2 <- RunUMAP(obj2, reduction = "pca", dims = 1:30)


saveRDS(obj2,file = "st.obj.rds")