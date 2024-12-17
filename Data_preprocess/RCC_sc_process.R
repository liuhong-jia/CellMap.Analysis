
sc.data <- readRDS("../data/RCC/TICAtlas_RNA.rds")
sc.obj <- subset(sc.data,subtype == "RCC")


sc.obj <-sc.obj %>% 
         NormalizeData() %>%
		 FindVariableFeatures() %>%
		 ScaleData() %>% 
		 RunPCA() %>% 
		 RunUMAP(dims = 1:30) %>%
		 FindNeighbors(., dims = 1:10)%>%
		 FindClusters(.,resolution = 0.8)

saveRDS(sc.obj,file = "../data/RCC/sc.obj.rds")