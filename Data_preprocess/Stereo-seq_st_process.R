
library(Seurat)
library(rhdf5)
library(Matrix)
cord_name='spatial'
meta <- read.table('../data/stereo/metadata.tsv',sep="\t",header=T,row.names=1)
pos <- read.table(paste0('../data/stereo/position_',cord_name,'.tsv'),sep="\t",header=T,row.names=1)


obj <- readRDS("../data/stereo/E16.5_E1S3_cell_bin_whole_brain.rds")

counts <- GetAssayData(obj,slot = "counts")
metadata <- meta
coords <- pos

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


st.obj <- createSpObj(counts, coord.df = coords, coord.label = c("x", "y"), meta.data = metadata)



##Extract subset object based on coordinates

st.sub <- subset(st.obj, image_x < -9600 &image_y < -2000, invert = FALSE)
coords <- st.sub@images$image@coordinates
Idents(st.sub) <- st.sub$region
st.obj <- subset(st.sub, subset = region %in% c("Pall", "SPall", "ChP","Pall VZ"))

st.obj <- SCTransform(st.obj, assay = "Spatial", ncells = 3000, verbose = FALSE)
st.obj <- RunPCA(st.obj)
st.obj <- RunUMAP(st.obj, dims = 1:30)
st.obj <- FindNeighbors(st.obj, dims = 1:30)
st.obj <- FindClusters(st.obj, resolution = 0.3, verbose = FALSE)

saveRDS(st.obj,file = "../data/stereo/st.obj.rds")
