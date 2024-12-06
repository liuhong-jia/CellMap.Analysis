##Estimate the global cell type proportions in ST data to generate CytoSPACE input file

cd /mnt/disk1/hongjia/work/data/VisiumHD/Seurat_v4/

library(Seurat)
library(dplyr)

sc.obj <- readRDS("sc.obj.rds")
st.obj <- readRDS("st.obj.rds")


ovp .genes <- intersect(rownames(sc.obj), rownames(st.obj))
sc.obj <- sc.obj[ovp.genes, ] %>% SCTransform(., assay = "RNA", verbose = TRUE)
sc.obj <- FindVariableFeatures(sc.obj)
st.obj <- st.obj[ovp.genes, ] %>% SCTransform(., assay = "Spatial", verbose = TRUE)
anchors <- FindTransferAnchors(reference = sc.obj, query = st.obj)

##celltype refers to the cell type lables column in scRNA-seq data 
predictions.assay <- TransferData(anchorset = anchors, refdata = sc.obj@meta.data[,"celltype"])
#st.obj <- AddMetaData(object = st.obj, metadata = predictions)

colnames(predictions.assay) <- gsub('prediction.score.', '', colnames(predictions.assay))
predictions.assay <- predictions.assay[, 2:(ncol(predictions.assay)-1)]
col_sums <- colSums(predictions.assay)
cellfrac <- col_sums / sum(col_sums)
cellfrac <- data.frame(Index = names(cellfrac), Fraction = cellfrac)


celltype <- names(table(sc.obj$celltype))

library(dplyr)
replace_celltype <- function(index_name, celltypes) {
  matched_type <- grep(index_name, celltypes, value = TRUE)
  if (length(matched_type) > 0) {
    return(matched_type[1])
  } else {
    return(index_name)  
  }
}

cellfrac$Index <- sapply(cellfrac$Index, replace_celltype, celltypes = celltype)

write.table(t(cellfrac),file = "your_path/Seurat_weights.txt",quote = FALSE, sep = "\t", col.names = FALSE)

