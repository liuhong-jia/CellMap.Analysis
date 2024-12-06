library(Seurat)
library(dplyr)

sc.data <- readRDS("your_path/sc.obj.rds")
st.data <- readRDS("your_path/st.obj.rds")

valid_types <- names(table(Idents(sc.data)))[table(Idents(sc.data)) >= 100]
sc.data_filtered <- subset(sc.data, idents = valid_types)

common_cells <- intersect(rownames(sc.data), rownames(st.data))
sc.data <- sc.data[common_cells, ]
st.data <- st.data[common_cells, ]

searchMarkersBySeurat <- function(obj, ...) {
  future::plan("multisession", workers = 9)
  deg.markers <- FindAllMarkers(obj, only.pos = TRUE, ...)
  return(split(deg.markers$gene, deg.markers$cluster))
}

markers <- searchMarkersBySeurat(sc.data)
top_100_markers <- lapply(markers, head, 100)
st.data <- AddModuleScore(st.data, features = top_100_markers, name = "CellTypeScore")
cell_type_names <- names(top_100_markers)
colnames(st.data@meta.data)[grep("CellTypeScore", colnames(st.data@meta.data))] <- cell_type_names
true.prop <- st.data@meta.data[, cell_type_names]


##Spot-level-correlation
  			 
common_columns <- intersect(colnames(true.prop), colnames(pre.prop))
true.prop <- true.prop[, common_columns, drop = FALSE]
pre.prop <- pre.prop[, common_columns, drop = FALSE]

true.prop <- true.prop[rownames(pre.prop), , drop = FALSE]

## Spot-level correlation

corre <- numeric(nrow(pre.prop))
for (i in 1:nrow(pre.prop)) {
  corre[i] <- cor(as.numeric(pre.prop[i, ]), as.numeric(true.prop[i, ]), use = "complete.obs")
}

corre_df <- data.frame(spot = rownames(pre.prop), correlation = corre)

## Normalized correlation
corre_df$corre <- scales::rescale(corre_df$corre, to = c(0, 1))			

	  