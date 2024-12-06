
##CellMap
library(dplyr)

##Read in the RCTD results and the mapping results of each tool.

RCTD <- readRDS("RCTD.rds")
CellMap <- readRDS("CellMap.rds")
metadata <- CellMap@meta.data[,c("CellType")]

colnames(RCTD) <- gsub("L2_3 IT", "L2/3 IT", colnames(RCTD))
prop <- prop.table(table(metadata[metadata %in% colnames(RCTD)]))
prop.RCTD <- colSums(RCTD[, names(prop)]) / sum(RCTD)
CellMap.R <- cor.test(prop, prop.RCTD, method = "pearson")




library(dplyr)

fraction <- read.csv("fractional_abundances_by_spot.csv")
RCTD <- readRDS("RCTD.rds")
rownames(fraction) <- fraction[, 1]
fraction <- fraction[, -1]

colnames(fraction) <- gsub("\\.", " ", colnames(fraction))
colnames(fraction) <- gsub("L2 3.IT", "L2/3 IT", colnames(fraction))
colnames(RCTD) <- gsub("L2_3 IT", "L2/3 IT", colnames(RCTD))
CellType <- c("Astro", "L2/3 IT", "L4", "L5 IT", "L5 PT", "NP", "L6 CT", "L6 IT", "L6b", "Pvalb", "Sst", "Vip", "Lamp5")
common_CellType <- intersect(CellType, intersect(colnames(RCTD), colnames(fraction)))
fraction <- fraction[, common_CellType]
RCTD <- RCTD[, common_CellType]
prop.CytoSpace <- colSums(fraction) / sum(fraction)
prop.RCTD <- colSums(RCTD) / sum(RCTD)
prop.RCTD <- prop.RCTD[names(prop.CytoSpace)]
CytoSPACE.R <- cor.test(prop.CytoSpace,prop.RCTD, method = "pearson")




###CellTrek
CellTrek <- readRDS("CellTrek.rds")
RCTD <- readRDS("RCTD.rds")

metadata <- CellTrek@meta.data[,c("cell_type")]
colnames(RCTD) <- sub("L2_3 IT", "L2/3 IT",colnames(RCTD))
CellType <- intersect(unique(metadata),colnames(RCTD))
metadata.1 <- metadata[metadata %in% CellType]

prop <- prop.table(table(metadata.1))

####RCTD
prop.RCTD <- colSums(RCTD[, CellType]) / sum(RCTD[, CellType])
prop.RCTD <- prop.RCTD[names(prop)]

CellTrek.R <- cor.test(prop.RCTD,prop, method = "pearson")


###Tangram

Tangram <- readRDS("Tangram.rds")
RCTD <- readRDS("RCTD.rds")
metadata <- Tangram@meta.data[,c("CellType")]

colnames(RCTD) <- sub("L2_3 IT", "L2/3 IT",colnames(RCTD))
CellType <- intersect(unique(metadata),colnames(RCTD))
metadata.1 <- metadata[metadata %in% CellType]

prop <- prop.table(table(metadata.1))

####RCTD
RCTD <- RCTD[,CellType]
prop.RCTD <- colSums(RCTD) / sum(RCTD)
prop.RCTD <- prop.RCTD[names(prop)]

Tangram.R <- cor.test(prop.RCTD,prop, method = "pearson")
