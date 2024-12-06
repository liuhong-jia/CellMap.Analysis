##MERFISH 

setwd("data/MERFISH/GSE113576/")
library(Seurat)
library(dplyr)
library(openxlsx)

obj.data <- Read10X(data.dir = "./")
obj <- CreateSeuratObject(obj.data,min.cells = 3,min.features = 200,project = "mouse_hypothalamus")
sc.obj <- obj %>% 
         NormalizeData() %>%
		     FindVariableFeatures() %>%
		     ScaleData() %>% 
		     RunPCA() %>% 
		     RunUMAP(dims = 1:30) %>%
		     FindNeighbors(., dims = 1:10)%>%
		     FindClusters(.,resolution = 0.8)

metadata <- read.xlsx("data/MERFISH/NIHMS1024025-supplement-Table_S1.xlsx",sheet = 1,startRow = 2)
metadata <- metadata[,c(1,2,3,4)]

rownames(metadata) <- metadata$Cell.name
metadata$Cell.name <- NULL
intersect(rownames(metadata), colnames(sc.obj))  
sc.obj@meta.data <- merge(sc.obj@meta.data, metadata, by = "row.names", all.x = TRUE)
rownames(sc.obj@meta.data) <- sc.obj@meta.data$Row.names
sc.obj@meta.data$Row.names <- NULL

head(sc.obj@meta.data)
colnames(sc.obj@meta.data)[8] <- "celltype"


saveRDS(sc.obj, file = "data/MERFISH/sc.obj.rds")


##Female
sc.obj.Female <- subset(sc.obj, subset = Sex == "Female")
Idents(sc.obj.Female) <- sc.obj.Female$celltype

sc.obj.Female <- subset(sc.obj.Female,idents = c("Ambiguous","Astrocytes","Endothelial",
                                     "Ependymal","Excitatory",
									 "Immature oligodendrocyte","Inhibitory",
									 "Microglia","Mature oligodendrocyte",
									 "Mural"))
									 
sc.obj.Female$celltype <- Idents(sc.obj.Female) %>% as.vector

sc.obj.Female$celltype[sc.obj.Female$celltype == "Astrocytes"] <- "Astrocyte"
sc.obj.Female$celltype[sc.obj.Female$celltype == "Immature oligodendrocyte"] <- "OD Immature"
sc.obj.Female$celltype[sc.obj.Female$celltype == "Mature oligodendrocyte"] <- "OD Mature"
sc.obj.Female$celltype[sc.obj.Female$celltype == "Mural"] <- "Mural/Pericyte"

table(sc.obj.Female$celltype)
Idents(sc.obj.Female) <- sc.obj.Female$celltype
sc.obj.Female <- subset(sc.obj.Female, celltype != "Ependymal")

downSamplSeurat <- function(obj, cluster.col = NULL, cnt = 200, seed = 123, percent = 0.3) {
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


sc.Female <- downSamplSeurat(sc.obj.Female,cnt = 2000,seed = 123,percent = NULL)

saveRDS(sc.Female,file = "data/MERFISH/sc.obj.Female.rds")

		
##Male
library(Seurat)
library(dplyr)

sc.obj.Male <- subset(sc.obj, subset = Sex == "Male")
Idents(sc.obj.Male) <- sc.obj.Male$celltype

sc.obj.Male <- subset(sc.obj.Male,idents = c("Ambiguous","Astrocytes","Endothelial",
                                     "Ependymal","Excitatory",
									 "Immature oligodendrocyte","Inhibitory",
									 "Microglia","Mature oligodendrocyte",
									 "Mural"))
									 
sc.obj.Male$celltype <- Idents(sc.obj.Male) %>% as.vector

sc.obj.Male$celltype[sc.obj.Male$celltype == "Astrocytes"] <- "Astrocyte"
sc.obj.Male$celltype[sc.obj.Male$celltype == "Immature oligodendrocyte"] <- "OD Immature"
sc.obj.Male$celltype[sc.obj.Male$celltype == "Mature oligodendrocyte"] <- "OD Mature"
sc.obj.Male$celltype[sc.obj.Male$celltype == "Mural"] <- "Mural/Pericyte"

table(sc.obj.Male$celltype)
Idents(sc.obj.Male) <- sc.obj.Male$celltype
sc.obj.Male <- subset(sc.obj.Male, celltype != "Ependymal")
sc.obj.Male <- subset(sc.obj.Male, celltype != "Ambiguous")

sc.Male <- downSamplSeurat(sc.obj.Male,cnt = 2000,seed = 123,percent = NULL)

saveRDS(sc.Male,file = "data/MERFISH/sc.obj.Male.rds")

		
		



