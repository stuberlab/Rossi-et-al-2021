rm(list = ls())  # Clear the environment
options(warn=-1) # Turn off warning message globally
library("Biobase")
library("Seurat")
library("ggplot2")
library("Matrix")
library(data.table)

objectDir = c("/Users/rossiadmin/Dropbox (Stuber Lab)/Mark/LHA projection paper/Data/Fig3 HCR")
getwd()
setwd(objectDir)
getwd()

genes = c("Vglut2","Vgat","Tomato","YFP","Nptx2","Pdyn","Pax6","Sostdc1","Pitx2")

## Big df with all HCR data ##
HCR.df = read.csv(paste0(objectDir, "/HCR_data.csv"))
HCR.df = as.data.frame(HCR.df, row.names=c(HCR.df$ID))
HCR.df

## make df with relevant metadata info #
metadata = as.data.frame(HCR.df[c('Sheet','Cell','Hemisphere','Mouse','Slice','MouseSlideSlice','Mouse.Hemisphere','Mouse.Hemisphere_Total_Cells','Fluor.Tag')],
                         row.names=c(HCR.df$ID))
tags=as.data.frame(metadata['Fluor.Tag'],colnames='Tag')
HCR.data <- HCR.df[genes]
HCR<- CreateSeuratObject(counts = t(HCR.data), min.cells = 0, min.features = 0, project = "LHA_HCR")
HCR <- AddMetaData(HCR, tags[1], col.name='Tag')


table(HCR@meta.data$Tag)


yfp.subset <- subset(x = HCR, subset = Tag == "YFP")
yfpcells <- WhichCells(yfp.subset)
tomato.subset <- subset(x = HCR, subset = Tag == "Tomato")
tomatocells <- WhichCells(tomato.subset)
fluor.subset <- subset(x = HCR, subset = (Tag == "Tomato" | Tag == "YFP"))
fluor.subset <- SetIdent(object=fluor.subset, cells=yfpcells,drop=TRUE,value='YFP' )
fluor.subset <- SetIdent(object=fluor.subset, cells=tomatocells,drop=TRUE, value='TOM' )
fluor.subset <- SubsetData(fluor.subset, ident.use=c("YFP","TOM"))
fluor.subset$ident <- Idents(fluor.subset)
fluor.subset

fluor.subset<- NormalizeData(object = fluor.subset,verbose = FALSE)
log1p<-Matrix(log(as.matrix(fluor.subset)+1), sparse = TRUE)


fluor.subset@assays$RNA@data<-t(log1p) # we do not do scaling so replace it simple log1p. important
fluor.subset <- ScaleData(fluor.subset,features=c("Vglut2","Vgat","Nptx2","Pdyn","Pax6","Sostdc1","Pitx2"), model.use='negbinom')

fluor.subset@assays$RNA@counts

# write all DEGs to csv
yfptomdeg_all_HCR <- FindMarkers(fluor.subset, ident.1 = "YFP", ident.2 = "TOM",min.pct=0,logfc.threshold=0,test.use='LR',slot="counts")
yfptomdeg_all_HCR
fwrite(x = yfptomdeg_all_HCR, row.names=TRUE, file = "HRC_seurat_deg_output.csv")
