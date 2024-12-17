##stereo-seq scRNA-seq preprocess

library(Seurat)
library(data.table)
library(Matrix)
mat <- readMM("../data/stereo/gene_sorted-matrix.mtx") 

mat <- as(mat,"dgCMatrix")
feature.names <- read.delim("../data/stereo/genes.tsv",header = FALSE,stringsAsFactors = FALSE)
barcode.names <- read.delim("../data/stereo/barcodes.tsv", header = FALSE,stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2

mat.1 <- mat[!duplicated(rownames(mat)), ]
obj <- CreateSeuratObject(counts = mat.1,project = "mouse.brain",assay = "RNA")
metadata <- read.table("../data/stereo/metaData_scDevSC.txt",sep = "\t",head = T,quote = "")
metadata <- metadata[-1,]
rownames(metadata) <- metadata$NAME
metadata <- metadata[,-1]
metadata <- metadata[, !(colnames(metadata) %in% c("Doublet_intersect", "Gral_cellType"))]

obj <- AddMetaData(obj, metadata = metadata)
sub.obj <- subset(obj,orig_ident %in% c("E16"))
Idents(sub.obj) <- sub.obj$New_cellType
counts <- table(sub.obj$New_cellType)
types <- names(counts[counts < 50])

sub.obj <- subset(sub.obj, idents = setdiff(Idents(sub.obj), types))
sub.obj$celltype <- sub.obj$New_cellType

cell_types<- c("Apical progenitors", "DL CPN", "Intermediate progenitors", 
                           "Interneurons", "Migrating neurons", "CThPN", "SCPN")

subset_obj <- subset(sub.obj, idents = cell_types)
subset_obj$celltype <- subset_obj$New_cellType

subset_obj <- subset_obj %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)

saveRDS(subset_obj,file = "../data/stereo/sc.obj.rds")