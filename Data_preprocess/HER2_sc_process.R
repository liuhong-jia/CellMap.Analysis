setwd("../data/GSE176078_Wu_etal_2021_BRCA_scRNASeq/Wu_etal_2021_BRCA_scRNASeq")

library(Seurat)
library(Matrix)
library(dplyr)

barcodes <- read.delim("count_matrix_barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes <- read.delim("count_matrix_genes.tsv", header = FALSE, stringsAsFactors = FALSE)
sparse_matrix <- readMM("count_matrix_sparse.mtx")

colnames(sparse_matrix) <- barcodes$V1
rownames(sparse_matrix) <- genes$V1

obj <- CreateSeuratObject(counts = sparse_matrix, project = "BRCA", min.cells = 3, min.features = 200)
metadata <- read.csv("metadata.csv")
rownames(metadata) <- metadata$X

metadata <- subset(metadata, select = -X)
obj <- AddMetaData(obj,metadata = metadata)

sc.obj <- subset(obj,subtype== "HER2+")
Idents(sc.obj) <- sc.obj$celltype_major

sc.obj <- sc.obj %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)
saveRDS(sc.obj, file = "../data/HER2/sc.obj.rds")

