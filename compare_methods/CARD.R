library(Seurat)
library(dplyr)
library(CARD)

sc.obj <- readRDS("your_path/brain_st_cortex.rds")
st.obj <- readRDS("your_path/brain_sc.rds")


st.counts <- GetAssayData(st.obj, slot = "count")
sc.counts <- GetAssayData(sc.obj, slot = "count")
sc.meta <- sc.obj@meta.data

## Set the cell type label column to cellType
sc.meta$cellType <- as.vector(sc.obj@meta.data[,"subclass"])
	
st.loc <- GetTissueCoordinates(object = obj.sp@images[[1]]) %>% `colnames<-`(c("x", "y"))
st.counts <- sweep(st.counts, 2, max(colSums(st.counts)), "/") * 1e6
st.counts <- apply(st.counts, 1, as.integer) %>% t()
colnames(st.counts) <- colnames(st.obj)

CARD.obj <- createCARDObject(
      sc_count = sc.counts,
      sc_meta = sc.meta,
      spatial_count = st.counts,
      spatial_location = st.loc,
      ct.varname = "cellType",
      ct.select = unique(sc.meta$cellType),
      sample.varname = "orig.ident",
      minCountGene = 100,
      minCountSpot = 5
    )
	
CARD.dec <- CARD_deconvolution(CARD_object = CARD.obj)

proportion = CARD.dec@Proportion_CARD

saveRDS(obj.sp,file = "your_path/CARD_decon.rds")
