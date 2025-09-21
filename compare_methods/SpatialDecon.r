library(Seurat)
library(SpatialDecon)
library(dplyr)

scrna_path <- args[1]
spatial_path <- args[2]
celltype_final <- args[3]
output_path <- args[4]


sc <- readRDS(scrna_path)
st <- readRDS(spatial_path)
  
colnames(sc@meta.data)[colnames(sc@meta.data) == celltype_final] <- "LabeledCellType"
sc$CellID <- rownames(sc@meta.data)
Idents(sc) <- sc$LabeledCellType
  
annots <- data.frame(cellType = as.character(Idents(sc)), 
                       cellID = names(Idents(sc)),
                       stringsAsFactors = FALSE)
  
# 构建signature matrix
custom_mtx_seurat <- create_profile_matrix(
    mtx = Seurat::GetAssayData(object = sc, assay = "RNA", slot = "counts"), 
    cellAnnots = annots, 
    cellTypeCol = "cellType", 
    cellNameCol = "cellID", 
    matrixName = paste0("custom_", dataset_name),
    outDir = NULL, 
    normalize = FALSE, 
    minCellNum = 5, 
    minGenes = 10
)

counts <- GetAssayData(st, slot= "counts")
andersson_g1 <- SeuratObject::CreateSeuratObject(counts = counts, assay="Spatial")

andersson_g1@meta.data$x <- st@meta.data$location_x
andersson_g1@meta.data$y <- st@meta.data$location_y
  
res <- runspatialdecon(andersson_g1, X = custom_mtx_seurat)
  
decon <- t(res$beta)
decon_prop <- decon / rowSums(decon)
  
write.csv(decon_prop, file = file.path(output_path, "spatialdecon_result.csv"))
