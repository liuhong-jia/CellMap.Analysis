
library(dplyr)

CellMap <- readRDS("CellMap.rds")

immune_cells <- filter(CellMap@meta.data, CellType %in% c("T/NK cell", "Myeloid cell","B cell"))

cell.counts <- immune_cells %>%
  group_by(SpotName) %>%
  summarise(ImmuneCellCount = n()) %>%
  as.data.frame


st.obj <- readRDS("st.obj.rds")
TLS_gene <- readxl::read_xlsx("TLS_gene.xlsx")

gene <- as.list(TLS_gene)
st.obj <- AddModuleScore(object = st.obj,
                         features = gene,
                         ctrl = 5,
                         name = "TLS_score"
                         )

score <- st.obj$TLS_score1 %>% as.vector 
st.obj$TLS <- scales::rescale(score, to = c(0, 1))			

TLS_score <- st.obj$TLS %>% as.data.frame
colnames(TLS_score) <- "score"
TLS_score$SpotName <- rownames(TLS_score)

data <- merge(TLS_score,cell.counts,by = "SpotName")
data$score_group <- cut(data$score, breaks = seq(0, 1, by = 0.2), include.lowest = TRUE)

corr_test <- cor.test(data$score,data$ImmuneCellCount, method = "pearson")