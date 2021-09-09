library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(matrixStats)
library(dplyr)
library(reshape2)
library(ape)
library(SeuratData)
library(cowplot)
library(patchwork)
library(gridExtra)
library(data.table)
library(magrittr)


################################################################################################################
################ Load data, select only neurons, set active assay to SCT and define features   #################

objectDir = c("/Users/rossiadmin/Desktop/seq")
setwd(objectDir)

# read objects
object <- readRDS(paste0(objectDir, "/200630_MLB017.subs.dbl.sct.int.pca.umap.neighbors40.minDist0.5.spread0.4.k15.res0.88.diet.neurons.integrated.rds"))

dim(object)

# subset to exclude cluster 1
object <- subset(object, idents = 1, invert=T)

DefaultAssay(object) <- "SCT"


## Features used to define clusters ##
features = c("Slc17a6","Slc32a1","tdtomato","eyfp","Nek7","Nptx2","Pcsk1","Rfx4","Hcrt","Pdyn","Lhx9","Plce1",
             "Pax6","Cartpt","Pmch","Aqp6","Ccdc42","Cyp19a1","Gsta2","Nppc","Sostdc1","Tcf4",
             "Nr4a2","Nfix","Prdm8","Onecut2","Onecut3","Trh","Ghrh","Ceacam10","Lamp5",
             "Qrfpr","Vipr2","St18","Crh","Sst","Gm47283","C1ql2","Otp","Lmo3","Lef1","Tcf7l2",
             "Foxd2","Ebf1","Foxb1","Nkx2???4","Cryab","Gpm6b","Sox2","Pgr15l","Meis2","Fezf2",
             "Scd2" ,"Ptn","Grpr","Lmo4","Fezf1","Tnnt2","Cabp7","Pitx2","Slc30a3")

## rename clusters based on above features ##
new.cluster.ids <- c("undef1","Hcrt","Glut1","Pmch1","Pmch2","GABA1","Glut2","Glut3","Glut4","Glut5",
                     "GABA2","GABA3","GABA4","GABA5","undef2","GABA6","Glut6","GABA7","GABA8","GABA9",
                     "Glut7","Glut8","Glut/GABA1","Glut9","GABA10","Glut10","Glut/GABA2","Glut/GABA3",
                     "GABA11","Glut11","Glut12","Glut13","Glut14","GABA12")
names(new.cluster.ids) <- levels(object)
object <- RenameIdents(object, new.cluster.ids)

features_hcr = c("Slc17a6","Slc32a1","tdtomato","eyfp","Nptx2","Pdyn","Pax6","Sostdc1","Pitx2")

features_subset = c("Slc17a6","Slc32a1","tdtomato","eyfp","Nptx2","Hcrt","Pdyn",
             "Pax6","Cartpt","Pmch","Gsta2","Nppc","Sostdc1","Tcf4",
             "Nfix","Trh","Lamp5","St18","Crh","Sst","Gm47283","C1ql2","Otp","Lmo3","Tcf7l2",
             "Foxd2","Foxb1","Gpm6b","Sox2","Pgr15l","Meis2",
             "Scd2" ,"Lmo4","Fezf1","Cabp7","Pitx2","Slc30a3")

################################################################################################################
##############################     Report median descriptive params                #############################

object@meta.data
UMIavg <- median(object@meta.data$nCount_SCT)
GENEavg <- median(object@meta.data$nFeature_SCT)
percentMitoavg <- median(object@meta.data$percentMito)
Readsavg <- median(object@meta.data$total_reads )
paste0("Number UMIs: ",UMIavg) 
paste0("Number GENES: ", GENEavg)
paste0("Percent Mito: ",percentMitoavg*100)
paste0("Number Reads: ",Readsavg)

################################################################################################################
############################### Dotplots for features of interest (all neurons)    #############################

savedir = c("/Users/rossiadmin/Dropbox (Stuber Lab)/Mark/LHA projection paper/Data/Fig2 seq")
setwd(savedir)

## Dotplots for feature ident ##
mindot=0.01
mincol=0
DotPlot(object, features = features,col.min=mincol,dot.min=mindot) + RotatedAxis()
ggsave2(paste0("Dotplot_all_dotmin",mindot,"_colmin",mincol,".pdf"), width=20,height=10,units="in",dpi=300)

DotPlot(object, features = features_hcr,col.min=mincol,dot.min=mindot) + RotatedAxis()
ggsave2(paste0("Dotplot_HCR_GENES_dotmin",mindot,"_colmin",mincol,".pdf"), width=20,height=10,units="in",dpi=300)

DotPlot(object, features = features_subset,col.min=mincol,dot.min=mindot) + RotatedAxis() +
  theme(axis.text.y = element_text(size = 18)) + theme(axis.text.x = element_text(size = 18))+
  theme(text=element_text(size=18))
ggsave2(paste0("Dotplot_subset_genes_dotmin",mindot,"_colmin",mincol,".pdf"), width=15,height=10,units="in",dpi=300)

################################################################################################################
########################### Load list of Tom/eYFP cells from Marcus ############################################

setwd(objectDir)
file = "210130_MLB017.subs.dbl.sct.int.pca.umap.neighbors40.minDist0.5.spread0.4.k15.res0.88.diet.neurons.integratedeyfp.tdtomato.cells.list.txt"
fluortable = read.table(file, header = FALSE, sep = "", dec = ".")
setwd(savedir)
yfptable = subset(fluortable, V4=='eyfp')
tomtable = subset(fluortable, V4=='tdtomato')
yfpandtomtable = subset(fluortable, V4=='eyfp/tdtomato')

yfpcells = yfptable[[1]]
tomcells = tomtable[[1]]
bothcells = yfpandtomtable[[1]]

dim(fluortable)
length(yfpcells)
length(tomcells)
length(bothcells)

################################################################################################################
###################### Parse YFP and TOM cells, plot, export data for python analysis ##########################

thresh=1

## save YFP and TOM data as csv ##
yfpdata <- subset(object, cells=yfpcells,features=features)
yfpdata <- as.matrix(GetAssayData(object = yfpdata, slot = "data",assay='SCT'))
yfpdata <-  t(yfpdata)
data_to_write_out <- as.data.frame(x = as.matrix(x = yfpdata))
fwrite(x = data_to_write_out, row.names=TRUE, file = "Seurat_SCT_YFP_Cells_subset_genes_UPDATED210130.csv")
# yfpdata <- GetAssayData(object, slot='data', assay='SCT')
# yfpdata <- t(object[features, yfpcells])
# write.table(yfpdata,'test.csv',sep = ',', row.names = F)


tomdata <- subset(object, cells=tomcells,features=features)
tomdata <- as.matrix(GetAssayData(object = tomdata, slot = "data",assay='SCT'))
tomdata <-  t(tomdata)
# tomdata <- GetAssayData(object, slot='data', assay='SCT')
# tomdata <- t(tomdata[features, tomcells])
data_to_write_out <- as.data.frame(x = as.matrix(x = tomdata))
fwrite(x = data_to_write_out, row.names=TRUE, file = "Seurat_SCT_Tomato_Cells_subset_genes_UPDATED210130.csv")


################################################################################################################
###################### Find % YFP and TOM cells expressing vglut2 or Vgat#######################################

### Vglut2 ######
vglut2cells <- WhichCells(object, expression = Slc17a6 > 1)
yfpandglutcells = intersect(yfpcells, vglut2cells)
tomandglutcells = intersect(tomcells, vglut2cells)

paste0("number of YFP+/Vglut2+ cells: ", length(yfpandglutcells))
paste0("Percent of YFP+/Vglut2+ cells: ", (length(yfpandglutcells)/length(yfpcells))*100)

paste0("number of Tom+/Vglut2+ cells: ", length(tomandglutcells))
paste0("Percent of Tom+/Vglut2+ cells: ", (length(tomandglutcells)/length(tomcells))*100)


### Vgat ######
vgatcells <- WhichCells(object, expression = Slc32a1 > 1)
yfpandvgatcells = intersect(yfpcells, vgatcells)
tomandvgatcells = intersect(tomcells, vgatcells)

paste0("number of YFP+/Vgat cells: ", length(yfpandvgatcells))
paste0("Percent of YFP+/Vgat cells: ", (length(yfpandvgatcells)/length(yfpcells))*100)

paste0("number of Tom+/Vgat cells: ", length(tomandvgatcells))
paste0("Percent of Tom+/Vgat cells: ", (length(tomandvgatcells)/length(tomcells))*100)

#### both ###
yfpvgatvglutcells= intersect(yfpandvgatcells, yfpandglutcells)
paste0("number of YFP+/Vgat+/vglut2+ cells: ", length(yfpvgatvglutcells))

tomvgatvglutcells= intersect(tomandvgatcells, tomandglutcells)
paste0("number of TOM+/Vgat+/vglut2+ cells: ", length(tomvgatvglutcells))

#### lepr #####
leprcells <- WhichCells(object, expression = Lepr > 0)
yfpandleprcells = intersect(yfpcells, leprcells)
tomandleprcells = intersect(tomcells, leprcells)
paste0("number of lepr+ cells: ", length(leprcells))
paste0("number of YFP+/lepr+ cells: ", length(yfpandleprcells))
paste0("number of YFP+/lepr- cells: ", length(yfpcells)-length(yfpandleprcells))
paste0("number of TOM+/lepr+ cells: ", length(tomandleprcells))
paste0("number of TOM+/lepr- cells: ", length(tomcells)-length(tomandleprcells))

#### Ghsr #####
ghsrcells <- WhichCells(object, expression = Ghsr > 0)
yfpandghsrcells = intersect(yfpcells, ghsrcells)
tomandghsrcells = intersect(tomcells, ghsrcells)
paste0("number of ghsr+ cells: ", length(ghsrcells))
paste0("number of YFP+/ghsr+ cells: ", length(yfpandghsrcells))
paste0("number of YFP+/ghsr- cells: ", length(yfpcells)-length(yfpandghsrcells))
paste0("number of TOM+/ghsr+ cells: ", length(tomandghsrcells))
paste0("number of TOM+/ghsr- cells: ", length(tomcells)-length(tomandghsrcells))


################################################################################################################
########################  get cluster assignments for YFP and Tom data   #######################################

yfpdata <- subset(object, cells=yfpcells,features=features)
yfpclustassignments = table(Idents(yfpdata))
write.table(x = yfpclustassignments, row.names=TRUE, sep=',', file = "YFP cluster assignments.csv")

tomdata <- subset(object, cells=tomcells,features=features)
tomclustassignments = table(Idents(tomdata))
write.table(x = tomclustassignments, row.names=TRUE, sep=',', file = "TOM cluster assignments.csv")


################################################################################################################
################ Compare YFP and TOM cells regardless of which neuronal cluster they come from #################

## create new seurat object comprising eYFP and Tomato expressing cells
# yfpandtom <- subset(object, subset=(eyfp > 0 | tdtomato > 0)) ##deprecated## 

yfpandtom <- subset(object, cells=c(yfpcells, tomcells,bothcells))
yfpandtom <- SetIdent(object=yfpandtom, cells=yfpcells,drop=TRUE,value='YFP' )
yfpandtom <- SetIdent(object=yfpandtom, cells=tomcells,drop=TRUE, value='TOM' )
yfpandtom <- SubsetData(yfpandtom, ident.use=c("YFP","TOM"))
yfpandtom$ident <- Idents(yfpandtom)
DefaultAssay(yfpandtom) <- "SCT"

yfpandtom

# write all DEGs to csv
yfptomdeg_all <- FindMarkers(yfpandtom, ident.1 = "YFP", ident.2 = "TOM",min.pct=.05,logfc.threshold=0,test.use='LR')
fwrite(x = yfptomdeg_all, row.names=TRUE, file = "YFPvTOM_DEGs_ALL_UPDATED210130.csv")

# subset degs based on FC and diff % for use here
yfptomdeg <- FindMarkers(yfpandtom, ident.1 = "YFP", ident.2 = "TOM",logfc.threshold = log(2),min.diff.pct=0.2)
degs<-rownames(yfptomdeg)

##save YFP and TOM data for plotting elsewhere
yfptomdeg_all_names <- rownames(yfptomdeg_all)
yfpandtomdata <- GetAssayData(yfpandtom, slot='data', assay='SCT') ## use SCT to generate output for SVM classifier
yfpandtomdata <- t(yfpandtomdata[yfptomdeg_all_names,])
data_to_write_out <- as.data.frame(x = as.matrix(x = yfpandtomdata))
fwrite(x = data_to_write_out, row.names=TRUE, file = "Seurat_YFP_and_Tomato_data_all_DEGs_UPDATED210430_test.csv")

RidgePlot(yfpandtom, features = degs, same.y.lims=T)
VlnPlot(yfpandtom,features=degs,split.by="ident",split.plot=T)


################################################################################################################
################################# save featureplots to multipage pdf     #######################################

figs <- list()
for(i in 1:length(features)){
  figs[[i]] <- FeaturePlot(object, features = features[i]) + xlim (c(-12, 3))
}
pdf('Featureplots.pdf')
for (p in figs){
  plot(p)
}
dev.off()


# # featureplots with no axes or legends
p <- FeaturePlot(object, features = features, blend = F, combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + theme(plot.title = element_text(size = 5)) + NoLegend() + NoAxes()+ xlim (c(-12, 3))
}
cowplot::plot_grid(plotlist = p)


################################################################################################################
################################# Make Feature plots for big cluster markers  ##################################

setwd(objectDir)
object <- readRDS(paste0(objectDir, "/200619_MLB017.subs.dbl.sct.int.pca.umap.neighbors40.minDist0.5.spread0.4.k15.res0.044.diet.sct.counts.rds"))
setwd(savedir)

cluster_markers <- c('Dcn','Igfbp7','Cx3cr1', 'Pdgfra', 'Gja1', 'Mag', 'Slc17a6', 'Slc32a1', 'Stmn2')

figs <- list()
for(i in 1:length(cluster_markers)){
  figs[[i]] <- FeaturePlot(object, features = cluster_markers[i], min.cutoff='q2')
}
pdf('Big_clusters_marker_genes_cutoff_q2.pdf')
for (p in figs){
  plot(p)
}
dev.off()

for(i in 1:length(cluster_markers)){
  plot<-FeaturePlot(object, features = cluster_markers[i], min.cutoff='q2', pt.size=3) + NoAxes()
  AugmentPlot(plot = plot,dpi=1000)
  ggsave(paste(cluster_markers[i],'.png'),plot=last_plot(),width=5,height=5, dpi=500, units="in")
}
