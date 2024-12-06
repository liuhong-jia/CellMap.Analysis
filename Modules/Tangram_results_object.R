# Generate the object for Tangram results.

library(openxlsx)

###Read in the Tangram running results:tangram.csv

metadata <- read.csv("tangram.csv")
metadata <- metadata[,-1]

##Read in the scRNA-seq and ST data

st.data <- readRDS("st.obj.rds")
sc.data <- readRDS("sc.obj.rds")
sc.data <- sc.data %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)

##Spatial coordinates of the ST slice
coord <- st.data@images$image@coordinates[,c(4,5)]
coord$SpotID <- rownames(coord)

data <- merge(coord,metadata,by = "SpotID")
data$UniqueCID <- seq(1:dim(data)[1])
locations <- data[,c(6,5,4,1,2,3)]
colnames(locations) <- c("UniqueCID","OriginalCID","CellType","SpotID","imagerow","imagecol")

library(dplyr)
num <- table(locations$SpotID)
unique_locations <- locations %>%
  select(SpotID, imagerow, imagecol) %>%
  distinct()
rownames(unique_locations) <- unique_locations$SpotID
spot.coords <- unique_locations[,-1]

getRandomCoords <- function(spot.coords,num) {
  #Estimated spot diameter
  kNN.dist <- dbscan::kNN(spot.coords,k= 4)$dist
  spot.diameter <- median(kNN.dist) %>% round
  spot.coords <- as.matrix(spot.coords)
  sc.coords <- lapply(1:nrow(spot.coords), function(i){
    sp.coord <- spot.coords[i,]
    cells <- num[i] %>% as.vector
    theta <- runif(cells,0,2*pi) 
    dis <- sqrt(runif(cells)) * (spot.diameter/2)
    cell.x <- sp.coord[1] + dis * cos(theta)
    cell.y <- sp.coord[2] + dis * sin(theta)
    coords <- data.frame(cell.x,cell.y)
    coords$SPOT <- paste0(spot.coords[i, 1], "x", spot.coords[i, 2])
    coords$centerX <- spot.coords[i, 1]
    coords$centerY <- spot.coords[i, 2]
    coords$SpotName <- rownames(spot.coords)[i]
    coords
  }) %>% do.call(rbind, .) %>% as.data.frame
  return (sc.coords)
}

coord.df <- getRandomCoords(spot.coords,num)
colnames(coord.df)[c(1,2,6)] <- c("imagerow","imagecol","SpotID") 

locations$imagerow <- coord.df$imagerow
locations$imagecol <- coord.df$imagecol

colnames(locations)[5] <- "cell.x"
colnames(locations)[6] <- "cell.y"

library(magrittr)
createSeuratObj <- function(st.data, sc.data,locations){
  
  sc.data$cell <- names(sc.data$orig.ident)
  
  obj <- CreateSeuratObject(counts = sc.data@assays$RNA[,locations$OriginalCID]%>% set_colnames(locations$UniqueCID) ,
                            project = 'tangram', assay = "RNA",
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
                                       assay = "RNA", key = 'tangram_')
  
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

saveRDS(obj,file = "your_path/Tangram.rds")