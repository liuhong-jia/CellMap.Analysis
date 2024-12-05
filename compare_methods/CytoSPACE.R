

##For low-resolution ST data

source('your_path/generate_cytospace_from_seurat_object.R')
scRNA_Seurat_Object <- readRDS("your_path/sc.obj.rds")
generate_cytospace_from_scRNA_seurat_object(scRNA_Seurat_Object, dir_out='.', fout_prefix='',rna_assay='RNA')
ST_Seurat_Object <- readRDS("your_path/st.obj.rds")
generate_cytospace_from_ST_seurat_object(ST_Seurat_Object, dir_out='.', fout_prefix='', slice='image')

cytospace -sp scRNA_data.txt -ctp cell_type_labels.txt -stp ST_data.txt -cp Coordinates.txt -o cytospace_results -sm lap_CSPR



##For high-resolution ST data
cytospace -sc \
   -sp scRNA_data.txt \
   -ctp cell_type_labels.txt \
   -stp ST_data.txt \
   -cp Coordinates.txt \
   -ctfep Seurat_weights.txt \
   -ncpsp cell.num.txt \
   -nop 2 