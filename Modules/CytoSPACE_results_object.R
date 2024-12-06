
library(Seurat)
library(dplyr)

##locations are the CytoSPACE results, containing the spatial coordinates of single cell.

locations <- read.csv("locations.csv")
locations <- locations[,c("UniqueCID","cell.x","cell.y")]


library(Seurat)
library(dplyr)

locations <- read.csv("locations.csv")

st.data <- readRDS("st.obj.rds")
sc.data <- readRDS("sc.obj.rds")

sc.data <- sc.data %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)
	
library(magrittr)


createSeuratObj <- function(st.data, sc.data,locations){
  
  sc.data$cell <- names(sc.data$orig.ident)
  
  obj <- CreateSeuratObject(counts = sc.data@assays$RNA[,locations$OriginalCID]%>% set_colnames(locations$UniqueCID) ,
                            project = 'CytoSPACE', assay = "RNA",
                            meta.data = sc.data@meta.data[locations$OriginalCID, ] %>%
                              mutate(Cell.new = locations$UniqueCID) %>%
                              set_rownames(locations$UniqueCID))
 
  colnames(locations)[1] <- "Cell.new"
  obj@meta.data <- dplyr::left_join(obj@meta.data,locations) %>% data.frame %>% set_rownames(obj$Cell.new)
 
  Idents(obj) <- obj$CellType
  
  obj[["RNA"]]@data <- obj[["RNA"]]@counts										
  sc.coord.obj <- CreateDimReducObject(embeddings = locations %>%
                                         dplyr::mutate(coord1 = cell.y, coord2 = max(cell.x)+ min(cell.x)- cell.x) %>%
                                         dplyr::select(c(coord1, coord2)) %>% set_rownames(locations$Cell.new) %>% as.matrix,
                                       assay = "RNA", key = 'CytoSPACE_')
  
  sc.pca.obj <- CreateDimReducObject(embeddings = sc.data@reductions$pca@cell.embeddings[locations$OriginalCID, ] %>%
                                       set_rownames(locations$Cell.new) %>% as.matrix, assay = "RNA", key = 'pca_')
  
  sc.umap.obj <- CreateDimReducObject(embeddings = sc.data@reductions$umap@cell.embeddings[locations$OriginalCID, ] %>%
                                        set_rownames(locations$Cell.new) %>% as.matrix, assay = "RNA", key = 'umap_')
  
  obj@reductions$CytoSpace <- sc.coord.obj
  obj@reductions$pca <- sc.pca.obj
  obj@reductions$umap <- sc.umap.obj
  
  if (!is.null(st.data)) {
    obj@images <- st.data@images
    obj@images[[1]]@assay <- DefaultAssay(obj)
    obj@images[[1]]@coordinates <- data.frame(imagerow = locations$cell.x, imagecol = locations$cell.y) %>% set_rownames(locations$Cell.new)
    obj@images[[1]]@scale.factors <- st.data@images[[1]]@scale.factors
  }
  return(obj)
}


obj <- createSeuratObj(st.data,sc.data,locations)

saveRDS(obj,file = "your_path/CytoSPACE.rds")