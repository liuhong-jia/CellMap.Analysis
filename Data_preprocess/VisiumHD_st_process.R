###Visium HD st data preprocess
library(Seurat)
library(hdf5r)

st.data <- Load10X_Spatial(data.dir = getwd(), bin.size = 8)
DefaultAssay(st.data) <- "Spatial.008um"

metadata <- read.csv("../data/VisiumHD/DeconvolutionResults_P1CRC.csv")

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


library(dplyr)

obj2 <- SCTransform(st.obj,assay = "Spatial") %>%
            RunPCA() %>%
            RunUMAP(dims = 1:30)%>%
			FindNeighbors()%>%
			FindClusters(resolution = 0.3)
           

saveRDS(obj2,file = "../data/VisiumHD/st.obj.rds")