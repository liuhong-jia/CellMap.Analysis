library(Matrix)
library(data.table)
library(Seurat)
library(dplyr)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(SeuratDisk)
library(future)
library(scran)

plan(multicore, workers = 8)
options(future.globals.maxSize = 3 * 1024 * 1024^2)


sc.data <- readRDS("your_path/sc.obj.rds")
st.data <- readRDS("your_path/st.obj.rds")


sce <- as.SingleCellExperiment(sc.data)
sce <- logNormCounts(sce)

spe <- as.SingleCellExperiment(st.data)
spe <- logNormCounts(spe)


genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(sce))
dec <- modelGeneVar(sce , subset.row = genes)

hvg <- getTopHVGs(dec, n = 3000)

colLabels(sce) <- colData(sce)$celltype

# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)

mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  x <- x[x$mean.AUC > 0.8, ]
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})

mgs_df <- do.call(rbind, mgs_fil)



### Deconvolution
res <- SPOTlight(
  x = counts(sce), 
  y = counts(spe),
  groups = as.character(sce$celltype_major), # 也可以是cluster，
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")
  
str(res) 
decon_mtrx <- res$mat

output_path <- "your_path"

write.csv(decon_mtrx, paste0(output_path, "/SPOTlight_decon.txt"))
