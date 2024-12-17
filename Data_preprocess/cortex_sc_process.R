###cortex

sc.obj <- readRDS("../data/cortex/allen_cortex.rds")

library(dplyr)
library(Seurat)

sc.obj <- sc.obj %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)

saveRDS(sc.obj,file = "../data/cortex/sc.obj.rds")
