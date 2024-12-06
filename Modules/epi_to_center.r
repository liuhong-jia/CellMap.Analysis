setwd("/mnt/disk1/hongjia/work/kidney/Epithelial_status/")

results <- readRDS("CellMap.rds")

## The central region object of the spot.
center.obj <- readRDS("center.rds")
center_data <- center.obj@meta.data[,c("imagerow","imagecol")]


#Euclidean distance
euclidean.dist <- function(vec1, vec2) {
  distance <- sqrt(sum((vec1 - vec2)^2))
  return(distance)
}


## The average Euclidean distance of cell types to the central region
calculate_average_distance <- function(cell_type_data, center_data) {
  num.data <- nrow(cell_type_data)
  num.center <- nrow(center_data)
  
  dist.matrix <- do.call(cbind, parallel::mclapply(1:num.data, function(i) {
    epi.vec <- cell_type_data[i,]
    dis <- sapply(1:num.center, function(j) {
      center.vec <- center_data[j,]
      euclidean.dist(epi.vec,center.vec)
    })
    dis
  }, mc.cores = parallel::detectCores()  - 1))
  
  average_distance <- mean(dist.matrix)
  
  return(average_distance)
}

cell_types <- unique(results$CellType) %>% as.vector

average_distances <- data.frame(CellType = character(0), AverageDistance = numeric(0))
for (cell_type in cell_types) {
  current_data <- subset(results, CellType == cell_type)
  current_data_coords <- current_data@meta.data[, c("cell.x", "cell.y")]
  average_distance <- calculate_average_distance(current_data_coords, center_data)
  average_distances <- rbind(average_distances, data.frame(CellType = cell_type, AverageDistance = average_distance))
}


## Data normalization
average_distances$AverageDistance <- scales::rescale(average_distances$AverageDistance, to = c(1, 50))

library(readxl)
truth <- read_xlsx("truth.xlsx") %>% as.data.frame()
colnames(truth) <- c("CellType","distance")

merged_data <- merge(average_distances, truth, by = "CellType")


r <- cor.test(merged_data$AverageDistance, merged_data$distance, method = "pearson")

