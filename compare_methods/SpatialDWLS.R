library(Seurat)  
library(ggplot2)
library(cowplot)
library(dplyr)
library(Giotto)
library(patchwork)
library(tidyverse)
library(Giotto)
library(SeuratDisk)

library(magick)
library(reticulate)

sc <- readRDS("your_path/sc.obj.rds")
st <- readRDS("your_path/st.obj.rds")

st_data <- createGiottoObject(
    raw_exprs = st@assays$Spatial@counts,
    #instructions = instrs
)

st_data <- normalizeGiotto(gobject = st_data)
st_data <- calculateHVG(gobject = st_data)

gene_metadata <- fDataDT(st_data)
featgenes <- gene_metadata[hvg == "yes"]$gene_ID
st_data <- runPCA(gobject = st_data, genes_to_use = featgenes, scale_unit = F)
signPCA(st_data, genes_to_use = featgenes, scale_unit = F)

st_data <- runUMAP(st_data, dimensions_to_use = 1:10)
st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 15)
st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000)

sc_data <- createGiottoObject(
    raw_exprs = sc@assays$RNA@counts,
    #instructions = instrs
)

sc_data <- normalizeGiotto(gobject = sc_data)
sc_data <- calculateHVG(gobject = sc_data)

gene_metadata <- fDataDT(sc_data)
featgenes <- gene_metadata[hvg == "yes"]$gene_ID
sc_data <- runPCA(gobject = sc_data, genes_to_use = featgenes, scale_unit = F)
signPCA(sc_data, genes_to_use = featgenes, scale_unit = F)
sc_data@cell_metadata$leiden_clus <- as.character(sc@meta.data[, "celltype"])
scran_markers_subclusters <- findMarkers_one_vs_all(
    gobject = sc_data,
    method = "scran",
    expression_values = "normalized",
    cluster_column = "leiden_clus"
)

Sig_scran <- unique(scran_markers_subclusters$genes[which(scran_markers_subclusters$ranking <= 100)])

norm_exp <- 2^(sc_data@norm_expr) - 1
id <- sc_data@cell_metadata$leiden_clus
ExprSubset <- norm_exp[Sig_scran, ]
Sig_exp <- NULL
for (i in unique(id)) {
    Sig_exp <- cbind(Sig_exp, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
}

colnames(Sig_exp) <- unique(id)
st_data <- runDWLSDeconv(st_data, sign_matrix = Sig_exp, n_cell = 5)

output_path = "your_path"

write.csv(st_data@spatial_enrichment$DWLS, paste0(output_path, "/SpatialDWLS_decon.txt"))