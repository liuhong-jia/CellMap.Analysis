###CellTrek
library(CellTrek)

sc.data <- readRDS("your_path/sc.obj.rds")
st.data <- readRDS("your_path/st.obj.rds")


####CellTrek
## Rename the cells/spots with syntactically valid names
obj <- RenameCells(st.data,new.names = make.names(Cells(st.data)))
sc.data <- RenameCells(sc.data, new.names=make.names(Cells(sc.data)))


traint <- CellTrek::traint(st_data= obj ,sc_data=sc.data, sc_assay='RNA', cell_names='celltype')

celltrek <- CellTrek::celltrek(st_sc_int=traint, int_assay='traint', sc_data = sc.data, sc_assay = 'RNA', 
                                   reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                   dist_thresh = 0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek

save(celltrek,file = "your_path/celltrek.rds")
