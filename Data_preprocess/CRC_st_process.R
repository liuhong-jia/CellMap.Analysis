
path <- "../data/CRC/st"

library(Matrix)
mat <- readMM(paste0(path,"GSM7058756_C1.matrix.mtx")) 
mat <- as(mat,"dgCMatrix")

feature.names <- read.delim(paste0(path,"GSM7058756_C1.features.tsv"),header = FALSE,stringsAsFactors = FALSE)
barcode.names <- read.delim(paste0(path,"GSM7058756_C1.barcodes.tsv"), header = FALSE,stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2

mat.1 <- mat[!duplicated(rownames(mat)), ]
obj <- CreateSeuratObject(counts = mat.1,project = "CRC",assay = "Spatial")


img <- Read10X_Image(image.dir = paste0(path,"spatial"))
DefaultAssay(object = img) <- "Spatial"
img <- img[colnames(x = obj)]
obj[["image"]] <- img

###标准化,降维和聚类
obj <- SCTransform(obj, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)


p <- SpatialDimPlot(obj, label = TRUE, label.size = 3)

saveRDS(obj,file = "../data/st.obj.rds")