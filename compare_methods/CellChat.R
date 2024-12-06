## Cell-cell interaction in Cellchat

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = FALSE)


##CellMap mapping results
results <- readRDS("CellMap.rds")

##Spatial transcriptomics data of renal cell carcinoma with tertiary lymphoid structures
st.obj <- readRDS("st.obj.rds")

##Cell-cell interactions in clusters with TLS structures

cluster_2 <- subset(st.obj,seurat_clusters == "2")
sub <- subset(results,SpotName %in% colnames(cluster_2))

cellchat <- createCellChat(sub@assays$RNA@data, meta = sub@meta.data, group.by = "CellType")
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object 

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multicore", workers = 20) # do parallel  
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat,PPI.human)  

# cellchat <- computeCommunProb(cellchat)
mycomputeCommunProb <- edit(computeCommunProb)  
environment(mycomputeCommunProb) <- environment(computeCommunProb)
cellchat <- mycomputeCommunProb(cellchat) 
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


pdf("inter.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("inter.1.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

 
pdf("network_plots.pdf", width = 10, height = 10)

mat <- cellchat@net$count
par(mfrow =c(3,3),xpd=T)

for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat),title.name = rownames(mat)[i])
}

dev.off()


## Visualize interactions with B cells separately
pdf("B_plots.pdf", width = 10, height = 10)
mat <- cellchat@net$count
mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
mat2[1,] <- mat[1,]
netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat),title.name = rownames(mat)[i])
dev.off()




##Cell-cell interactions in clusters without TLS structures

cluster_0 <- subset(st.obj,seurat_clusters == "0")
sub.1 <- subset(results,SpotName %in% colnames(cluster_0))


cellchat <- createCellChat(sub.1@assays$RNA@data, meta = sub.1@meta.data, group.by = "CellType")
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object 

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multicore", workers = 20) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  

# cellchat <- computeCommunProb(cellchat) 
mycomputeCommunProb <- edit(computeCommunProb) 
environment(mycomputeCommunProb) <- environment(computeCommunProb)
cellchat <- mycomputeCommunProb(cellchat) 
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

pdf("inter.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("inter.1.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("network_plots.pdf", width = 10, height = 10)

mat <- cellchat@net$count
par(mfrow =c(3,3),xpd=T)

for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat),title.name = rownames(mat)[i])
}

dev.off()



## Visualize interactions with B cells separately

pdf("B_plots.pdf", width = 10, height = 10)
mat <- cellchat@net$count
mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
mat2[1,] <- mat[1,]
netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat),title.name = rownames(mat)[i])
dev.off()