library(Seurat)
library(dplyr)
library(SeuratDisk)
library(future)

input_path <- "../data/"
datasets <- list.dirs(input_path, full.names = TRUE, recursive = FALSE)

for (dataset_dir in datasets) {
  dataset_name <- basename(dataset_dir)
  cat("Processing:", dataset_name, "\n")
  result_file <- file.path(dataset_dir, "true.matrix.rds")
  if (file.exists(result_file)) {
    cat("Result already exists for", dataset_name, "- skipping\n")
    next
  }

  sc_file <- file.path(dataset_dir, "scRNA.h5seurat")
  st_file <- file.path(dataset_dir, "Spatial.h5seurat")
  
  if (!file.exists(sc_file) || !file.exists(st_file)) {
    cat("Missing scRNA.h5seurat or Spatial.h5seurat in", dataset_name, "\n")
    next
  }

  sc.data <- LoadH5Seurat(sc_file)
  st.data <- LoadH5Seurat(st_file)

  if (!"celltype_final" %in% colnames(sc.data@meta.data)) {
    cat("'celltype_final' not found in", dataset_name, "- skipping\n")
    next
  }
  
  Idents(sc.data) <- sc.data$celltype_final
  cell_counts <- table(Idents(sc.data))
  valid_types <- names(cell_counts[cell_counts >= 50])
  sc.data <- subset(sc.data, idents = valid_types)

  inter <- intersect(rownames(sc.data), rownames(st.data))
  sc.data <- sc.data[inter, ]
  st.data <- st.data[inter, ]

  searchMarkersBySeurat <- function(obj, ...) {
    options(future.globals.maxSize = 40000 * 1024^2)
    future::plan("multisession", workers = 10)
    deg.markers <- FindAllMarkers(obj, only.pos = TRUE,...)	
    return(split(deg.markers$gene, deg.markers$cluster))
  }

  markers <- searchMarkersBySeurat(sc.data)
  top_100_markers <- lapply(markers, function(cell_type) head(cell_type, 100))

  st.data <- AddModuleScore(object = st.data, features = top_100_markers, name = "CellTypeScore",ctrl = 80,nbin = 20)

  cell_type_names <- names(top_100_markers)
  names(st.data@meta.data)[grep("CellTypeScore", names(st.data@meta.data))] <- cell_type_names
  true.matrix <- st.data@meta.data[, cell_type_names]

  saveRDS(true.matrix, file = result_file)
  cat("Saved:", result_file, "\n")
}

cat(">>> All datasets processed.\n")