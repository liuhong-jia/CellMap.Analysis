library(spacexr)
library(Seurat)

reference <- readRDS("your_path/sc.obj.rds")

###celltype cell type column in scRNA-seq data.
Idents(reference) <- reference$celltype
counts <- table(Idents(reference))
filtered <- names(counts[counts >= 100])
reference <- subset(reference, idents = filtered)

counts <- GetAssayData(reference,slot = "counts") 
meta_data <- reference@meta.data 
cell_types <- meta_data$celltype; names(cell_types) <- rownames(meta_data) 
cell_types <- as.factor(cell_types)

reference <- Reference(counts, cell_types)
print(dim(reference@counts)) 
table(reference@cell_types) 
saveRDS(reference, file = "your_path/SCRef.rds")


st.data <- readRDS("st.obj.rds")
counts <- GetAssayData(st.data,slot= "counts")
coords <- st.data@images$image@coordinates[,c(4,5)]

puck <- SpatialRNA(coords, counts)
saveRDS(puck,file = "your_path/puck.rds")
print(dim(puck@counts)) 

myRCTD <- create.RCTD(puck, reference, max_cores = 10)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
str(myRCTD)

results <- myRCTD@results

norm_weights = normalize_weights(results$weights)

saveRDS(norm_weights,file = "your_path/RCTD_decon.rds")
