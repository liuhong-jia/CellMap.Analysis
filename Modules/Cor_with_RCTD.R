
##CellMap
setwd("../work/cortex/")

library(dplyr)
RCTD <- readRDS("RCTD/RCTD.rds")
sc.out <- readRDS("CellMap/sc.out.rds")
sub <- subset(sc.out,idents = c("Astro","L2/3 IT","L4","L5 IT","L5 PT","NP","L6 CT","L6 IT","L6b","Pvalb","Sst","Vip","Lamp5","Sncg"))
metadata <- sub@meta.data[,c("CellType")]

colnames(RCTD) <- gsub("L2_3 IT", "L2/3 IT", colnames(RCTD))
prop <- prop.table(table(metadata[metadata %in% colnames(RCTD)]))
prop.RCTD <- colSums(RCTD[, names(prop)]) / sum(RCTD)
CellMap.R <- cor.test(prop, prop.RCTD, method = "pearson")



library(dplyr)
fraction <- read.csv("cytospace_results/fractional_abundances_by_spot.csv")
RCTD <- readRDS("RCTD/RCTD.rds")
rownames(fraction) <- fraction[, 1]
fraction <- fraction[, -1]

colnames(fraction) <- gsub("\\.", " ", colnames(fraction))
colnames(fraction) <- gsub("L2 3.IT", "L2/3 IT", colnames(fraction))
colnames(RCTD) <- gsub("L2_3 IT", "L2/3 IT", colnames(RCTD))

CellType <- c("Astro","L2/3 IT","L4","L5 IT","L5 PT","NP","L6 CT","L6 IT","L6b","Pvalb","Sst","Vip","Lamp5","Sncg")
fraction <- fraction[,CellType]
RCTD <- RCTD[,CellType]

prop.CytoSpace <- colSums(fraction) / sum(fraction)
prop.RCTD <- colSums(RCTD) / sum(RCTD)
prop.RCTD <- prop.RCTD[names(prop.CytoSpace)]
CytoSPACE.R <- cor.test(prop.CytoSpace,prop.RCTD, method = "pearson")


###CellTrek
CellTrek <- readRDS("CellTrek/results.rds")
RCTD <- readRDS("RCTD/RCTD.rds")

metadata <- CellTrek@meta.data[,c("cell_type")]
colnames(RCTD) <- sub("L2_3 IT", "L2/3 IT",colnames(RCTD))

CellType <- c("Astro","L2/3 IT","L4","L5 IT","L5 PT","NP","L6 CT","L6 IT","L6b","Pvalb","Sst","Vip","Lamp5","Sncg")
metadata.1 <- metadata[metadata %in% CellType]
prop <- prop.table(table(metadata.1))

####RCTD
prop.RCTD <- colSums(RCTD[, CellType]) / sum(RCTD[, CellType])
prop.RCTD <- prop.RCTD[names(prop)]

CellTrek.R <- cor.test(prop.RCTD,prop, method = "pearson")


###Tangram

Tangram <- readRDS("Tangram/obj.rds")
RCTD <- readRDS("RCTD/RCTD.rds")

CellType <- c("Astro","L2/3 IT","L4","L5 IT","L5 PT","NP","L6 CT","L6 IT","L6b","Pvalb","Sst","Vip","Lamp5","Sncg")
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
