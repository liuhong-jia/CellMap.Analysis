##########################################################################################################
#' @title downSamplSeurat
#' @description Sample single cells by cell type.
#' @param obj Seurat object.
#' @param cnt Sample size for each ident, default: 200.
#' @param seed Randon seed, default: 123.
#' @return Subset of seurat object.
#' @export

##########################################################################################################

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

##########################################################################################################
#' @title perturb_sc_seura
#' @description Random perturbation on single-cell transcriptomic data
#' @param sc.obj Seurat object of scRNA-seq data
#' @param ratio The proportion of gene perturbation,set to 0,0.05 and 0.2.

##########################################################################################################	

#' @example sc.train <- perturb_sc_seurat(sc.train, ratio = 0.1)


perturb_sc_seurat <- function(sc.obj, ratio = 0.1) {
   sc.matrix <- sc.obj$RNA@data %>% as.matrix
   num.genes <- nrow(sc.matrix)
   perturb.genes.count <- round(num.genes * ratio)
    perturb.genes <- sample(1:num.genes, perturb.genes.count)
   sc.matrix[perturb.genes, ] <- sc.matrix[perturb.genes, ][sample(length(perturb.genes)), ]
   sc.obj$RNA@data <- sc.matrix
   return(sc.obj)
}

##########################################################################################################
#' @title getInterdata
#' @description Integrate single-cell(sc) and spatial transcriptomic(st) data.
#' @param sc.data sc data Seurat object.
#' @param st.data st data Seurat object.
#' @param nfeatures number of feature genes for intergration.
#' @param coord coordinates column names in st images slot.
#' @return Seurat object of integrated sc and st data.
#' @export

##########################################################################################################	

#' @example sc.st.obj <- getInterdata(sc.train, st.data, nfeatures = 2000,coord = c("imagerow","imagecol"))

getInterdata <- function(sc.data, st.data, nfeatures = 2000,coord = c("imagerow","imagecol")) {

  sc.data$id <- colnames(sc.data)
  st.data$id <- colnames(st.data)
  
  sc.data$type <- "sc"
  st.data$type <- "st"
  
  st.data <- st.data[,rownames(st.data@images[[1]]@coordinates)]
  
  st.data$image_x <- st.data@images[[1]]@coordinates[, coord[1]]
  st.data$image_y <- st.data@images[[1]]@coordinates[, coord[2]]
  
  DefaultAssay(st.data) <- "Spatial"
  DefaultAssay(sc.data) <- "RNA"
  
  # Select integration features
  
  sc.st.features <- SelectIntegrationFeatures(list(st.data,sc.data), nfeatures = nfeatures)
  
  # Find transfer anchors
  sc.st.anchors <- FindTransferAnchors(reference = sc.data, query = st.data,
                                               features = sc.st.features, reduction = "cca")
  
  # Transfer data
  st.data.trans <- TransferData(anchorset = sc.st.anchors, refdata = GetAssayData(sc.data,slot = "data")[sc.st.features,],weight.reduction = "cca")
  
  st.data@assays$trans <- st.data.trans
  
  # Create integrated Seurat object
  sc.st.meta <- bind_rows(st.data@meta.data, sc.data@meta.data)
  counts <- cbind(data.frame(st.data@assays$trans@data), data.frame(sc.data@assays$RNA@data[sc.st.features,]))
  
  rownames(sc.st.meta) <- sc.st.meta$id
  colnames(counts) <- sc.st.meta$id
  
  sc.st.obj <- CreateSeuratObject(counts = counts, assay = "integrate", meta.data = sc.st.meta)
  
  sc.st.obj@assays$integrate@data <- sc.st.obj@assays$integrate@counts
  sc.st.obj@assays$integrate@counts <- matrix(nrow = 0,ncol = 0)
  
  # Run PCA and UMAP
  sc.st.obj <- ScaleData(sc.st.obj,features = rownames(sc.st.obj)) %>% RunPCA(features = rownames(sc.st.obj)) %>% RunUMAP(dim = 1:30)
  
  sc.st.obj@images <- st.data@images
  sc.st.obj@images[[1]]@coordinates <- data.frame(imagerow = sc.st.obj$image_x,
                                                  imagecol = sc.st.obj$image_y)
  return(sc.st.obj)
}

##########################################################################################################
#' @title getDistMatrix
#' @description Calculate the distance matrix between single cells and spatial spots.
#' @param sc.st.obj Seurat object of integrated sc and st data.
#' @return A distance matrix of single cell and spot.
#' @export
##########################################################################################################	
#' @example dist <- getDistMatrix(sc.st.obj)

getDistMatrix <- function(sc.st.obj) {
  
    UMAP <- sc.st.obj@reductions$umap@cell.embeddings
	
	st.meta.data <- sc.st.obj@meta.data[sc.st.obj@meta.data$type == "st",]
	sc.meta.data <- sc.st.obj@meta.data[sc.st.obj@meta.data$type == "sc",]
	
	st.UMAP <- UMAP[rownames(st.meta.data),]
	sc.UMAP <- UMAP[rownames(sc.meta.data),]
    
    #euclidean distance
	euclidean.dist <- function(vec1, vec2) {
	distance <- sqrt(sum((vec1 - vec2)^2))
	return(distance)
	}	
	num.cells <- nrow(sc.UMAP)
	num.spots <- nrow(st.UMAP)
	dist.matrix <- matrix(0, nrow = num.cells, ncol = num.spots)
	
	#Parallel computing
    num.cores <- parallel::detectCores()
	
	dist.matrix <- do.call(cbind, parallel::mclapply(1:num.cells, function(i) {
	cell.vec <- sc.UMAP[i,]
		dis <- sapply(1:num.spots, function(j) {
			spot.vec <- st.UMAP[j,]
			euclidean.dist(cell.vec, spot.vec)
		})
		dis
	}, mc.cores = num.cores - 1))
	
	rownames(dist.matrix) <- rownames(st.UMAP)
	colnames(dist.matrix) <- rownames(sc.UMAP)
	
	return(dist.matrix)
}

##########################################################################################################
#' @title getCells
#' @description Calculate the 5 nearest single cells for each spot based on the distance matrix between single cells and spots.
#' @param sc.st.dist A matrix representing the spatial distance between single cells and spots.
#' @param k The number of nearest neighboring single cells.
#' @return A list corresponding to the nearest neighboring single cells for each spot.
#' @export
##########################################################################################################	
#' @example nearCells <- getCells(dist, k = 5)


getCells <- function(sc.st.dist, k = 5) {
  neighbors <- lapply(1:nrow(sc.st.dist), function(i){
    near.index <- order(sc.st.dist[i,])[1:k]
    colnames(sc.st.dist)[near.index]
  })
  names(neighbors) <- rownames(sc.st.dist)
  return(neighbors)
}

##########################################################################
#' @title SynSpatialData 
#' @description Generate simulated spatial transcriptome data.
#' @param nearCells A list of the 5 nearest neighbors of each spot.
#' @param sc.data sc data Seurat object.
#' @param st.data st data Seurat object.
#' @param ratio Percentage of perturbed genes.Defalut:0.05.
#' @return A matrix of simulated spatial transcriptome data.
#' @export

##########################################################################
#' @example spot.data <- SynSpatialData(nearCells,sc.train,st.data,ratio = 0.05)


SynSpatialData <- function(nearCells, sc.data, st.data) {

  sc.matrix <- GetAssayData(sc.data, slot = "data") %>% as.matrix

  numSpots <- length(nearCells) 
  numGenes <- nrow(sc.matrix)
  
  combinedData <- matrix(0, nrow = numGenes, ncol = numSpots)
  colnames(combinedData) <- names(nearCells)

  # Define a function for parallel processing
  process_spot <- function(spotId) {
    neighborIds <- nearCells[[spotId]]
    tempData <- rowSums(sc.matrix[, neighborIds])
    return(tempData)
  }
  
  # Use mclapply for parallel processing
  nearCellNames <- names(nearCells)
  numCores <- detectCores()
  tempDataList <- parallel::mclapply(nearCellNames, process_spot, mc.cores = numCores - 2)
  
  # Combine the temporary data into combinedData
  for (i in 1:length(nearCells)) {
    combinedData[,i] <- tempDataList[[i]]
  }
  rownames(combinedData) <- rownames(sc.data)
  colnames(combinedData) <- names(nearCells)
  
  # Add perturbation
 
  return(combinedData)
}

##########################################################################
#' @title createSeuratObj 
#' @description Create Seurat object.
#' @param spot.data A matrix of simulated spot expression data.
#' @param st.data st data Seurat object.
#' @return A spatial transcriptome Seurat object.
#' @export

##########################################################################
#' @example obj <- createSeuratObj(spot.data,st.data)

createSeuratObj <- function(spot.data, st.data) {
  obj <- CreateSeuratObject(counts = spot.data, project = 'simulation', assay = "Spatial")
  obj@images <- st.data@images

  obj <- SCTransform(obj, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
  obj <- RunPCA(obj, assay = "SCT", verbose = FALSE,npcs = 50)
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:10)
  obj <- FindClusters(obj, resolution = 0.001, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:10)
  
  SpatialDimPlot(obj, label = TRUE, label.size = 3)
  return(obj)
}



library(Seurat)
library(dplyr)
library(parallel)
st.data <- readRDS("your_path/st.obj.rds")
reference <- readRDS("your_path/sc.obj.rds")

##Celltype refers to the cell type column in the single-cell reference.
Idents(reference) <- reference$celltype

## Split the single-cell reference into training and testing sets in a 1:1 ratio, with the training set used to generate simulated spatial transcriptomics data.

sc.train <- downSamplSeurat(reference,percent = 0.5)

##perform random perturbation on single-cell transcriptomic data.
sc.train <- perturb_sc_seurat(sc.train, ratio = 0.05)
test.cell <- setdiff(colnames(reference),colnames(sc.train))
sc.test <- reference[,test.cell]

sc.st.obj <- getInterdata(sc.train, st.data, nfeatures = 2000,coord = c("imagerow","imagecol"))
dist <- getDistMatrix(sc.st.obj)
nearCells <- getCells(dist, k = 5)
spot.data <- SynSpatialData(nearCells, sc.train, st.data)
simu.obj <- createSeuratObj(spot.data,st.data)


saveRDS(simu.obj ,file = "your_path/simu.obj.rds")
saveRDS(sc.test ,file = "your_path/ref.sc.obj.rds")
