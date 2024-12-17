##Processing of spatial transcriptomics data for the first slice

data.st <- read.csv("../data/MERFISH/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv")
data.st.bregma0.26 <- data.st[data.st$Bregma==0.26,]
table(data.st.bregma0.26$Animal_ID)
data.st.bregma0.26.id1 <- data.st.bregma0.26[data.st.bregma0.26$Animal_ID==1,]


count.data <- data.st.bregma0.26.id1[,10:170]
rownames(count.data) <- data.st.bregma0.26.id1$Cell_ID
count.data <- t(count.data)
count.data2 <- na.omit(count.data)

metadata <- data.st.bregma0.26.id1[,1:9]
rownames(metadata) <- data.st.bregma0.26.id1$Cell_ID

createSpObj <- function(counts, coord.df, coord.label = c("x", "y"), meta.data = metadata, class = "SlideSeq") {
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

df <- metadata[,c("Centroid_X","Centroid_Y")]
colnames(df) <- c("x","y")

st.obj <- createSpObj(counts = count.data2,coord.df= df,meta.data = metadata)

st.obj <- SCTransform(st.obj, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
st.obj <- FindVariableFeatures(st.obj, assay = "Spatial", selection.method = "vst", nfeatures = 2000)
st.obj <- ScaleData(st.obj,assay = "Spatial") 
st.obj <- RunPCA(st.obj,verbose = FALSE)
st.obj <- FindNeighbors(st.obj, reduction = "pca", dims = 1:30)
st.obj <- FindClusters(st.obj, verbose = FALSE,resolution = 0.5)
st.obj <- RunUMAP(st.obj, reduction = "pca", dims = 1:30)

st.obj@meta.data <- st.obj@meta.data %>%
  dplyr::select(-Cell_ID, -Centroid_X, -Centroid_Y, -Neuron_cluster_ID)

colnames(st.obj@meta.data)[8] <- "celltype"

##Merge cell subtypes into broader cell types
celltype_mapping <- c(
  "Endothelial 1" = "Endothelial",
  "Endothelial 2" = "Endothelial",
  "Endothelial 3" = "Endothelial",
  "OD Immature 1" = "OD Immature",
  "OD Immature 2" = "OD Immature",
  "OD Mature 1" = "OD Mature",
  "OD Mature 2" = "OD Mature",
  "OD Mature 3" = "OD Mature",
  "OD Mature 4" = "OD Mature"
)

library(dplyr)

st.obj$celltype <- st.obj$celltype %>%
  recode(!!!celltype_mapping)

st.obj$celltype[st.obj$celltype == "Pericytes"] <- "Mural/Pericyte"
st.obj <- subset(st.obj, celltype != "Ependymal")
st.obj <- subset(st.obj, celltype != "Ambiguous")


saveRDS(st.obj,file = "../data/MERFISH/st.obj.rds")



