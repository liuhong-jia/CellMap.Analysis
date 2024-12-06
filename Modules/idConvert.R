#' @title Gene id transformation of ST data.
#' @description Convert the ensembl gene id into gene symbol of the ST data.
#' @param scMap_obj An scMap object.
#' @return An scMap object.
#' @export
scMap.id.convert <- function(scMap_obj){
  ST.data <- scMap_obj@input$counts
  if (grep("ENSG",rownames(ST.data))){
    symbol <- clusterProfiler::bitr(rownames(ST.data), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
    ST.data <- ST.data[symbol[,1],]
    rownames(ST.data) <- symbol[,2]
    scMap_obj@input$counts <- ST.data
  }
  else {
    scMap_obj@input$counts <- ST.data
  }
  scMap_obj
}