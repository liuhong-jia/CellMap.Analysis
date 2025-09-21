library(Seurat)
library(dplyr)

args <- commandArgs(T)
scrna_path <- args[1]
spatial_path <- args[2]
celltype_final <- args[3]
output_path <- args[4]

sc_rna <- readRDS(scrna_path)
sc_rna <- SCTransform(sc_rna, verbose = FALSE)
spatial <- readRDS(spatial_path)
spatial <- SCTransform(spatial, verbose = FALSE)

anchors <- FindTransferAnchors(reference = sc_rna, query = spatial, dims = 1:30, normalization.method = "SCT")

if (!(celltype_final %in% colnames(sc_rna@meta.data))) {
    cat("Column", celltype_final, "not found in", dataset_name, "\n")
    next
}

predictions <- TransferData(anchorset = anchors, refdata = sc_rna@meta.data[,celltype_final], dims = 1:30)
write.csv(predictions, file = file.path(output_path, "Seurat_result.txt"))