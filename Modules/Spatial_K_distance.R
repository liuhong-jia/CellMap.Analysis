

##Computation of spatial k-distance

kdist <- function(inp_df, ref=NULL, ref_type='all', que=NULL, k = 10, new_name='kdist', keep_nn=F) {
  ## Check ##
  if (!all(c('CellType', 'coord_x', 'coord_y') %in% colnames(inp_df))) {stop("Input data must contain 'CellType', 'coord_x', 'coord_y' columns")}
  if (any(!c(ref %in% inp_df$CellType, que %in% inp_df$CellType))) {stop('Reference or query group not in CellType')}

  que_dat <- inp_df[inp_df$CellType %in% que, c('coord_x', 'coord_y')]
  kNN_dist_df <- data.frame(matrix(NA, nrow=nrow(que_dat), ncol=0))
  kNN_nn_list <- list()
  if (ref_type=='each' & length(ref)>0) {
    for (i in 1:length(ref)) {
      ref_i <- ref[i]
      ref_dat <- inp_df[inp_df$CellType %in% ref_i, c('coord_x', 'coord_y')]
      ## error when k >= nrow(ref_dat), reset k ##
      if (nrow(ref_dat) <= k) {k <- (nrow(ref_dat)-1)}

      kNN_res <- dbscan::kNN(x=ref_dat, k=k, query=que_dat)
      kNN_dist <- apply(kNN_res$dist, 1, mean)
      kNN_dist_df <- base::cbind(kNN_dist_df, kNN_dist)
      if (keep_nn) {
        nn_mat <- kNN_res$id
        nn_mat[] <- rownames(ref_dat)[c(nn_mat)]
        kNN_nn_list[[i]] <- nn_mat
      } else {
        kNN_nn_list[[i]] <- matrix(NA, 0, 0)
      }
    }
    colnames(kNN_dist_df) <- paste0(ref, '_kdist')
    names(kNN_nn_list) <- paste0(ref, '_ref')
  } else {
    ref_dat <- inp_df[inp_df$CellType %in% ref, c('coord_x', 'coord_y')]
    kNN_res <- dbscan::kNN(x=ref_dat, k=k, query=que_dat)
    kNN_dist_df <- data.frame(apply(kNN_res$dist, 1, mean)) %>% magrittr::set_colnames(new_name)
    if (keep_nn) {
      nn_mat <- kNN_res$id
      nn_mat[] <- rownames(ref_dat)[c(nn_mat)]
      kNN_nn_list$nn_ref <- nn_mat
    } else {
      kNN_nn_list[[1]] <- matrix(NA, 0, 0)
    }
  }
  output <- list(kdist_df=kNN_dist_df, knn_list=kNN_nn_list)
  return(output)
}

setwd("/mnt/disk1/hongjia/work/cortex/CellMap/")

CellMap <- readRDS("CellMap.rds")
sub <- subset(CellMap,idents = c("Astro","L2/3 IT","L4","L5 IT","L5 PT","NP","L6 IT","L6 CT","L6b"))
CellMap_data <- sub@meta.data[,c("CellType","cell.x","cell.y")]
colnames(CellMap_data) <- c("CellType","coord_x","coord_y")
rownames(CellMap_data) <- NULL
test_df <- CellMap_data

kdist_out_CellMap <- kdist(inp_df=test_df, ref="Astro", ref_type='each', que=unique(test_df$CellType), k=10, keep_nn=T)
CellMap_data <- cbind(CellMap_data,kdist_out_CellMap$kdist_df)




##CellTrek
CellTrek <- readRDS("CellTrek.rds")
Idents(CellTrek) <- CellTrek$cell_type
sub <- subset(CellTrek,idents = c("Astro","L2/3 IT","L4","L5 IT","L5 PT","NP","L6 IT","L6 CT","L6b"))
CellTrek_data <- sub@meta.data[,c("cell_type","coord_x","coord_y")]
colnames(CellTrek_data) <- c("CellType","coord_x","coord_y")
test_df <- CellTrek_data

kdist_out_CellTrek <- kdist(inp_df=test_df, ref="Astro", ref_type='each', que=unique(test_df$CellType), k=10, keep_nn=T)
CellTrek_data <- cbind(CellTrek_data,kdist_out_CellTrek$kdist_df)

saveRDS(CellTrek_data,file = "your_path/CellTrek_result.rds")


## CytoSPACE
setwd("/mnt/disk1/hongjia/work/cortex/data/cortex/cytospace_results/")
locations <- read.csv("locations.csv")
CytoSPACE_data <- locations[,c("CellType","cell.x","cell.y")]
cell_types <- c("Astro","L2/3 IT","L4","L5 IT","L5 PT","L6 IT","L6 CT","L6b")
CytoSPACE_data <- subset(CytoSPACE_data, CellType %in% cell_types)
colnames(CytoSPACE_data) <- c("CellType","coord_x","coord_y")

test_df <- CytoSPACE_data

kdist_out_CytoSPACE <- kdist(inp_df=test_df, ref="Astro", ref_type='each', que=unique(test_df$CellType), k=10, keep_nn=T)
CytoSPACE_data <- cbind(CytoSPACE_data,kdist_out_CytoSPACE$kdist_df)

saveRDS(CytoSPACE_data,file = "your_path/CytoSPACE_result.rds")


## Tangram
Tangram <- readRDS("Tangram.rds")

Idents(Tangram) <- Tangram$cell_type
sub <- subset(Tangram,idents = c("Astro","L2/3 IT","L4","L5 IT","L5 PT","NP","L6 IT","L6 CT","L6b"))
Tangram_data <- sub@meta.data[,c("cell_type","cell.x","cell.y")]
colnames(Tangram_data) <- c("CellType","coord_x","coord_y")

test_df <- Tangram_data
kdist_out_Tangram <- kdist(inp_df=test_df, ref="Astro", ref_type='each', que=unique(test_df$CellType), k=10, keep_nn=T)
Tangram_data <- cbind(Tangram_data,kdist_out_Tangram$kdist_df)

saveRDS(Tangram_data,file = "your_path/Tangram_result.rds")
