
rm(list=ls())
setwd("../data/CRC/")
library(data.table)

expr <- fread("GSM7058754_immune_counts.txt",sep = "\t",head = T)
metadata <- read.table("GSM7058754_immune_meta.txt",sep = "\t",head = T)
expr[1:3,1:3]
genes <- expr$V1
expr_matrix <- as.data.frame(expr)
rownames(expr_matrix) <- genes
expr_matrix[1:3,1:3]
expr_matrix <- expr_matrix[,-1]

library(Seurat)
obj <- CreateSeuratObject(counts = expr_matrix, project = "CRC", min.cells = 3, min.features = 200)
rownames(metadata) <- metadata$X
metadata <- subset(metadata, select = -X)
obj <- AddMetaData(obj,metadata = metadata)

library(dplyr)
immune.obj <- obj

##重新归类细胞类型
immune.obj$CellType <- NA

immune.obj$CellType[immune.obj$cluster == "B01_plasma_IGKC"] <- "B cell"
immune.obj$CellType[immune.obj$cluster == "B02_plasma_IGLC3"] <- "B cell"
immune.obj$CellType[immune.obj$cluster == "B03_B_MS4A1"] <- "B cell"
immune.obj$CellType[immune.obj$cluster == "B04_plasma_IGHG2"] <- "B cell"
immune.obj$CellType[immune.obj$cluster == "B05_plasma_IGLL1"] <- "B cell"
immune.obj$CellType[immune.obj$cluster == "B06_B_IGHD"] <- "B cell"
immune.obj$CellType[immune.obj$cluster == "B07_plasma_PPIB"] <- "B cell"
immune.obj$CellType[immune.obj$cluster == "B08_plasma_IGLJ3"] <- "B cell"
immune.obj$CellType[immune.obj$cluster == "B09_GCB_RGS13"] <- "B cell"

immune.obj$CellType[immune.obj$cluster == "M01_Mono_CD14"] <- "Myeloid cell"
immune.obj$CellType[immune.obj$cluster == "M02_Mac_CXCL9"] <- "Myeloid cell"
immune.obj$CellType[immune.obj$cluster == "M03_cDC_CD1C"] <- "Myeloid cell"
immune.obj$CellType[immune.obj$cluster == "M04_Mac_ZNF331"] <- "Myeloid cell"
immune.obj$CellType[immune.obj$cluster == "M05_Mac_SPP1"] <- "Myeloid cell"
immune.obj$CellType[immune.obj$cluster == "M06_Mac_HSPA6"] <- "Myeloid cell"
immune.obj$CellType[immune.obj$cluster == "M07_Mac_SERPINB2"] <- "Myeloid cell"
immune.obj$CellType[immune.obj$cluster == "M08_cycling_MKI67"] <- "Myeloid cell"
immune.obj$CellType[immune.obj$cluster == "M09_cDC_CPNE3"] <- "Myeloid cell"
immune.obj$CellType[immune.obj$cluster == "M10_Mono_FCGR3A"] <- "Myeloid cell"
immune.obj$CellType[immune.obj$cluster == "M11_cDC_LAMP3"] <- "Myeloid cell"
immune.obj$CellType[immune.obj$cluster == "Mast cell"] <- "Myeloid cell"

immune.obj$CellType[immune.obj$cluster == "N01_NK_NCAM1"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "N02_NK_FCGR3A"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "pDC"] <- "pDC"
immune.obj$CellType[immune.obj$cluster == "T01_cycling_TNK"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T02_gdT_TRGC1"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T03_CD8_CCL4"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T04_CD8_CD69"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T05_CD8_HSPA6"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T06_CD8_ITGA1"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T07_CD8_CCL20"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T08_CD8_CX3CR1"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T09_CD8_TXNIP"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T10_CD8_CXCL13"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T11_CD4_BHLHE40"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T12_CD4_SELL"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T13_CD4_HSPA6"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T14_Treg_FOXP3"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T15_CD4_ZBTB10"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T16_Treg_TNFRSF9"] <- "T/NK cell"
immune.obj$CellType[immune.obj$cluster == "T17_CD4_CXCL13"] <- "T/NK cell"

###非免疫细胞
library(data.table)
expr <- fread("GSM7058755_non_immune_counts.txt",sep = "\t",head = T)
metadata <- read.table("GSM7058755_non_immune_meta.txt",sep = "\t",head = T)

metadata$X <- gsub("-",".",metadata$X)


expr[1:3,1:3]
genes <- expr$V1
expr_matrix <- as.data.frame(expr)
rownames(expr_matrix) <- genes
expr_matrix[1:3,1:3]
expr_matrix <- expr_matrix[,-1]

library(Seurat)
obj <- CreateSeuratObject(counts = expr_matrix, project = "CRC", min.cells = 3, min.features = 200)
rownames(metadata) <- metadata$X
metadata <- subset(metadata, select = -X)
obj <- AddMetaData(obj,metadata = metadata)

library(dplyr)
no.immune.obj <- obj

##重新归类细胞类型
no.immune.obj$CellType <- NA

no.immune.obj$CellType[no.immune.obj$cluster == "E01_endothelial_SELP"] <- "Endothelial cell"
no.immune.obj$CellType[no.immune.obj$cluster == "E02_endothelial_DLL4"] <- "Endothelial cell"
no.immune.obj$CellType[no.immune.obj$cluster == "E03_endothelial_NOTCH3"] <- "Endothelial cell"
no.immune.obj$CellType[no.immune.obj$cluster == "E04_endothelial_CD36"] <- "Endothelial cell"
no.immune.obj$CellType[no.immune.obj$cluster == "E05_cycling_MKI67"] <- "Endothelial cell"
no.immune.obj$CellType[no.immune.obj$cluster == "E06_endothelial_CLEC4G"] <- "Endothelial cell"
no.immune.obj$CellType[no.immune.obj$cluster == "F01_fibroblast_PRELP"] <- "Fibroblast"
no.immune.obj$CellType[no.immune.obj$cluster == "F02_fibrblast_MCAM"] <- "Fibroblast"
no.immune.obj$CellType[no.immune.obj$cluster == "F03_fibroblast_CXCL14"] <- "Fibroblast"
no.immune.obj$CellType[no.immune.obj$cluster == "F04_fibroblast_C3"] <- "Fibroblast"
no.immune.obj$CellType[no.immune.obj$cluster == "F05_fibroblast_COCH"] <- "Fibroblast"
no.immune.obj$CellType[no.immune.obj$cluster == "F06_cycling_MKI67"] <- "Fibroblast"

no.immune.obj$CellType[no.immune.obj$cluster == "Tu01_AREG"] <- "Tumor cell"
no.immune.obj$CellType[no.immune.obj$cluster == "Tu02_DEFA5"] <- "Tumor cell"
no.immune.obj$CellType[no.immune.obj$cluster == "Tu03_SRRM2"] <- "Tumor cell"
no.immune.obj$CellType[no.immune.obj$cluster == "Tu04_RGMB"] <- "Tumor cell"
no.immune.obj$CellType[no.immune.obj$cluster == "Tu05_PCNA"] <- "Tumor cell"
no.immune.obj$CellType[no.immune.obj$cluster == "Tu06_NKD1"] <- "Tumor cell"
no.immune.obj$CellType[no.immune.obj$cluster == "Tu07_MKI67"] <- "Tumor cell"
no.immune.obj$CellType[no.immune.obj$cluster == "Tu08_GNG13"] <- "Tumor cell"
no.immune.obj$CellType[no.immune.obj$cluster == "Tu09_MUC2"] <- "Tumor cell"
no.immune.obj$CellType[no.immune.obj$cluster == "Tu10_COL3A1"] <- "Tumor cell"
no.immune.obj$CellType[no.immune.obj$cluster == "Tu11_PLA2G2A"] <- "Tumor cell"

sc.obj <- merge(immune.obj,no.immune.obj)
saveRDS(sc.obj,file = "/mnt/disk1/hongjia/work/CRC/data/sc.obj")


sc.obj$CellType[sc.obj$cluster == "pDC"] <- "pDC"


downSamplSeurat <- function(obj, cluster.col = NULL, cnt = 3000, seed = 123, percent = 0.1) {
  set.seed(seed)
  if (!is.null(cluster.col)) Idents(obj) <- obj@meta.data[, cluster.col]
  cells <- Idents(obj) %>% table
  sub.cells <- sapply(names(cells), function(xx) {
    sub.cells <- Idents(obj)[Idents(obj) == xx] %>% names
    cnt <- ifelse(is.null(percent), cnt, length(sub.cells) * percent)
    if (length(sub.cells) > cnt) sub.cells <- sample(sub.cells, cnt, replace = FALSE)
    return(sub.cells)
  }) %>% unlist(use.names = F)
  subset(obj, cells = sub.cells)
}

Idents(sc.obj) <- sc.obj$CellType

sc.obj.5 <- downSamplSeurat(sc.obj,cluster.col = "CellType",cnt = 5000, seed = 123,percent = NULL)
sc.obj.5 <- subset(sc.obj.5,idents = "pDC",invert = T)

sc.obj.5 <- sc.obj.5 %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)
saveRDS(sc.obj.5,file = "../data/CRC/sc.obj.obj")
