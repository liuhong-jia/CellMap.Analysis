#' @title simulate_spatial_data
#' @description Generate simulated spatial transcriptome data.
#' @param num_spots The number of spots generated. Default: 1000.
#' @param obj A scRNA-seq data object.
#' @param seed Randon seed, default: 1234.
#' @return The gene expression matrix of the spot and the cell type composition of each spot.
#' @export

#' @example simu.data <- simulate_spatial_data(num_spots = 1000,obj,seed = 1234)

simulate_spatial_data <- function(num_spots = 1000,obj,seed = 1234) { 
  cell_types <- levels(obj)
  expression_levels <- GetAssayData(obj,slot = "counts")
  cell_type_counts <- replicate(num_spots, sample(cell_types, size = sample(3:9),replace = FALSE))
  num_cells_per_spot <- replicate(num_spots, sample(10:50, size = 1))
  spots <- list()
  spot_cell_types <- list()
  for (i in 1:num_spots) {
    set.seed(seed)
    cells <- Idents(obj) %>% table
    expression <- matrix(0, nrow = nrow(expression_levels), ncol = 1)
    types <- c()
    for (j in 1:num_cells_per_spot[i]) {
      cell_type <- sample(cell_type_counts[[i]], size = 1)
      types <- c(types, cell_type)
      sub.cells <- sample(Idents(obj)[Idents(obj) == cell_type],size = 1) %>% names
      cell_expression <- expression_levels[,sub.cells] %>% as.matrix
      spot_expression <- expression + cell_expression
      #exp <- spot_expression
      #/num_cells_per_spot[i]
    }
    spots[[i]] <-  spot_expression
    spot_cell_types[[i]] <- types
  } 
  
  proportions <- lapply(spot_cell_types, function(x) {
    counts <- table(x)
    prop <- counts / sum(counts)
    prop[setdiff(all_cell_types, names(prop))] <- 0
    return(prop)
  })
  
  proportions_matrix <- do.call(rbind, proportions) %>% t %>% `colnames<-`(seq(1,num_spots))
  spot <- do.call(cbind,spots) %>% `colnames<-`(seq(1,num_spots))
  
  return(list(spots = spot,spot_cell_types = spot_cell_types,proportion = proportions_matrix))
}


#' @title synSpotData 
#' @description Generate simulated spatial transcriptome data.
#' @param sc.data sc data Seurat object.
#' @param st.data st data Seurat object.
#' @param lambda Expected value of Poisson distribution.
#' @param perturb.ratio Percentage of perturbed genes.Defalut:0.05.
#' @return A list A list containing the synthetic spot expression matrix, the number of single cells per spot, and the correspondence between spots and cells.
#' @export

#' @example spot.data <- synSpotData(sc.data,st.data,lambda = 5,ratio = 0.05)

synSpotData <- function(sc.data,st.data,lambda = 5,ratio = 0.05) {
	num.genes <- nrow(sc.data)
	num.iter = ncol(st.data)
	syn.data <- matrix(0,nrow = num.genes,ncol = num.iter)
    values <- numeric(num.iter)
	
	cell.names.list <- list()
	
	for (i in 1:num.iter) {
		random <- rpois(length(lambda),lambda)
		if(random == 0) {random =1}
		values[i] <- random
		matrix <- GetAssayData(sc.data,slot = "counts")
		cells <- sample(1:ncol(matrix), random, replace = TRUE)
		cell.names <- colnames(matrix[,cells]) 
		cell.names.list[[i]] <- cell.names
		
		if(values[i] <= 1){
		    data <- matrix[, cells] %>% as.matrix
		}else{
		    data <- rowSums(matrix[, cells]) %>% as.matrix
        }
        syn.data[, i] <- data 
	}
	rownames(syn.data) <- rownames(sc.data)
	perturb <- round(nrow(syn.data) * ratio)
	perturb.genes <- sample(1:nrow(syn.data),perturb)
	# Replace gene names of perturbed genes with randomly selected gene names
	syn.data[perturb.genes, ] <- as.matrix(syn.data[perturb.genes, ])[sample(length(perturb.genes)), ]
	colnames(syn.data) <- colnames(st.data)
	names(cell.names.list) <- colnames(syn.data)
	
	df.list <- lapply(names(cell.names.list), function(spot.name) {
    cell.names <- cell.names.list[[spot.name]]
    if (length(cell.names) > 0) {
        df <- data.frame(spot = spot.name, cell = cell.names)
        return(df)
    } else {
        return(NULL) 
    }
    })
    df.list <- df.list[!sapply(df.list, is.null)]
	
    final.df <- do.call(rbind, df.list)
	syn.data <- syn.data %>% data.frame
	return(list(syn.data = syn.data, num = values, spot.to.cells = final.df))
}







