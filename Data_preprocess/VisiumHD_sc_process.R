## Visium HD scRNA-seq preprocess

library(Seurat)
library(dplyr)
crc.obj <- Read10X_h5("HumanColonCancer_Flex_Multiplex_count_filtered_feature_bc_matrix.h5")


crc.sc.obj <- CreateSeuratObject(crc.obj)
crc.sc.obj$orig.ident <- "CRC"
metadata <- read.csv("SingleCell_MetaData.csv")
rownames(metadata) <- metadata$Barcode  
crc.sc.obj <- AddMetaData(crc.sc.obj, metadata = metadata)

p1.sc.obj <- subset(crc.sc.obj,subset = Patient=="P1CRC")
p1.sc.obj <- subset(p1.sc.obj,subset = QCFilter=="Keep")

sc.obj <- p1.sc.obj %>% 
         NormalizeData() %>%
		 FindVariableFeatures() %>%
		 ScaleData() %>% 
		 RunPCA() %>% 
		 RunUMAP(dims = 1:30) %>%
		 FindNeighbors(., dims = 1:10)%>%
		 FindClusters(.,resolution = 0.8)

Idents(sc.obj) <- sc.obj$Level1
sc.obj$celltype <- sc.obj$Level1



###Rename cell types

level2_to_level1 <- c(
  "CAF" = "Fibroblast",
  "CD4 T cell" = "T cells",
  "CD8 Cytotoxic T cell" = "T cells",
  "Endothelial" = "Endothelial",
  "Enteric Glial" = "Neuronal",
  "Enterocyte" = "Intestinal Epithelial",
  "Epithelial" = "Intestinal Epithelial",
  "Fibroblast" = "Fibroblast",
  "Goblet" = "Intestinal Epithelial",
  "Lymphatic Endothelial" = "Endothelial",
  "Macrophage" = "Myeloid",
  "Mast" = "Myeloid",
  "Mature B" = "B cells",
  "mRegDC" = "Myeloid",
  "Myofibroblast" = "Fibroblast",
  "Neuroendocrine" = "Neuronal",
  "Neutrophil" = "Myeloid",
  "pDC" = "Myeloid",
  "Pericytes" = "Endothelial",
  "Plasma" = "B cells",
  "Proliferating Immune II" = "T cells",
  "SM Stress Response" = "Smooth Muscle",
  "Smooth Muscle" = "Smooth Muscle",
  "Tuft" = "Intestinal Epithelial",
  "Tumor I" = "Tumor",
  "Tumor II" = "Tumor",
  "Tumor III" = "Tumor",
  "Tumor V" = "Tumor",
  "Unknown III (SM)" = "Smooth Muscle",
  "vSM" = "Smooth Muscle"
)

sc.obj$celltype <- sc.obj$Level2
sc.obj$celltype <- level2_to_level1[sc.obj$Level2]
Idents(sc.obj) <- sc.obj$celltype


saveRDS(sc.obj,file = "sc.obj.rds")

