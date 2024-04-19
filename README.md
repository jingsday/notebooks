# notebooks

###Read Cell Rangell output into objects

##R Seurat pipeline
library(Seurat)
#work directory and three files location
setwd("")
cells<- "filtered_barcodes.tsv.gz"
features<- "features.tsv.gz"
mtx <-"matrix.mtx.gz"#

#Read three files into one matrix
data <-ReadMtx(mtx=mtx, cells=cells,features=features,feature.column=2)

#Create a seurat object for further analysis
seurat_obj<-CreateSeuratObject(counts = data, project = "name", min.cells = 3)

##Python Scanpy pipeline
import pandas as pd
import scanpy as sc

#Read documents into adata 
mtx ="matrix.mtx.gz"
adata = sc.read_mtx(mtx)

cells=pd.read_csv('barcodes.tsv.gz',header=None)
features=pd.read_csv('features.tsv.gz',header=None,sep='\t')
adata= adata.T
adata.obs['CellID']= cells
adata.var['Gene']= features[1].tolist()
adata.var.index= adata.var['Gene']
adata.var_names_make_unique() 
