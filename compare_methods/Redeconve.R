library(Redeconve)

sc <- readRDS("your_path/sc.obj.rds")
st <- readRDS("your_path/st.obj.rds")

library(Seurat)
sc.data <- GetAssayData(sc,slot = "counts")
st.data <- GetAssayData(st,slot = "counts")

## celltype refers to the cell type labels in the scRNA-seq data

annotations <- sc@meta.data[,c("celltype")] %>% as.data.frame
colnames(annotations) <- "annotations"

annotations$barcodes <- rownames(sc@meta.data)
annotations <- annotations[,c(2,1)]

## Coordinates of the spatial transcriptomics slice
coords <- st@images$image@coordinates[,c(4,5)]

ref <- get.ref(sc.data, annotations, dopar = F)

## deconvolution
res.ct <- deconvoluting(ref, st.data, genemode = "def", hpmode = "auto", dopar = T, ncores = 10)


write.table(res.ct,file = "your_path/Redeconve.decon.txt",sep = "\t",quote = T)

proportions <- sweep(t(res.ct), 1, rowSums(t(res.ct)), FUN = "/")
write.table(proportions,file = "your_path/Redeconve.result.txt",sep = "\t",quote = T)
