# references:
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# https://satijalab.org/seurat/articles/integration_mapping.html
# https://satijalab.org/seurat/articles/integration_introduction.html
library(Seurat)
library(dplyr)
#library(magrittr)
#library(hdf5r)
library(ggplot2)


################################################
####  10X datasets pre-processing:  ####
################################################

# -- load 10X dataset
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P33")
P33.d <- Read10X(data.dir = "P33")

# -- Create seurat object (non-normalized data)
## --- features: detected in > 3 cells
## --- cells: > 200 features are detected
P33<-CreateSeuratObject(count= P33.d, project = "P33", min.cells = 3, min.features = 200)

# -- QC checking:
## --- calculate mitochondrial percent:
P33[["percent.mt"]] <- PercentageFeatureSet(P33, pattern = "^mt-")


## --- label samples by development time:
P33[["time"]]  <-   "P33"


## --- Visualize QC metrics as a violin plot:
VlnPlot(P33, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(P33, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(P33, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# -- QC, clean datasets: (NOTE: key step, change threshold accordingly)
Rmin <- 200
Rmax <- 5000
mt.per <-10
P33 <- subset(P33, subset = nFeature_RNA > Rmin & nFeature_RNA < Rmax & percent.mt < mt.per)
P33 # check cells


#################################################################
#### analysis for directly merged datasets: (for comparison) ####
#################################################################

s.object <- P33  # define Seurat object

# -- Normalizing the data
s.object <- NormalizeData(s.object, normalization.method = "LogNormalize", scale.factor = 10000)

# -- Identification of highly variable features 
s.object <- FindVariableFeatures(s.object, selection.method = "vst", nfeatures = 2000)#nfeatures = 2000:	Number of features to select as top variable features
top20<- head(VariableFeatures(s.object), 20)
top20

plot1 <- VariableFeaturePlot(s.object)
plot2 <- LabelPoints(plot = plot1, points = top20)
plot2
plot1 + plot2

# -- Scaling the data
all.genes <- rownames(s.object) # get gene names
s.object<- ScaleData(s.object, features = all.genes)

# -- Perform linear dimensional reduction
s.object <- RunPCA(s.object,npcs = 50, features = VariableFeatures(object = s.object)) # default: 50 PCs

## --- Examine and visualize PCA results a few different ways (optional)
print(s.object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(s.object, dims = 1:5, reduction = "pca")
DimHeatmap(s.object, dims = 1:15, cells = 500, balanced = TRUE) #cells = 500:seleceted 500 cells to run pca

## --- show PCA scores:  for determining the ??dimensionality??
npca <- 20  # define number of PCA scores
### opt-1:
ElbowPlot(s.object,ndims = npca)
### opt-2: NOTE: The following process can take a LONG ~~ time
s.object <- JackStraw(s.object,dims = npca, num.replicate = 100) #num.replicate:Number of replicate samplings to perform
s.object <- ScoreJackStraw(s.object, dims = 1:npca)
JackStrawPlot(s.object, dims = 1:npca)
#optional: output PCA done object 
#saveRDS(s.object, file = "merged_only_PCAdone.object")
#s.object <- readRDS("merged_only_PCAdone.object")

# -- Cluster the cells
ndim <-20  # dims according to 'significant' PCs
## --- Computes the k.param nearest neighbors for a given dataset
s.object <- FindNeighbors(s.object, dims = 1:ndim)
## --- change resolution to get different clusters. (between 0.4-1.2 typically)
s.object <- FindClusters(s.object, resolution =0.5) 
head(Idents(s.object), 5) # check cluster IDs of the first 5 cells

# -- Run non-linear dimensional reduction 
s.object <- RunUMAP(s.object, dims = 1:ndim)
p1 <- DimPlot(s.object, reduction = "umap",label = TRUE)
p1
p1_time <- DimPlot(s.object, reduction = "umap",label = TRUE,group.by = "time")
p1_time
p1+p1_time
saveRDS(s.object, file = "P33.object")

# -- find markers for every cluster compared to all remaining cells (postive & negtive)
s.object.markers <- FindAllMarkers(s.object , min.pct = 0.25, logfc.threshold = 0.25)
#memory.limit()
#memory.limit(1000000)
# -- plotting the top 20 markers (or all markers if less than 20) for each cluster
top20 <- s.object.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
DoHeatmap(s.object, features = top20$gene)+scale_fill_gradientn(colors = c("skyblue", "white", "#FF3300"))
# -- output makers
cluster.n <- length(levels(s.object.markers$cluster)) # get number of clusters
write.table(s.object.markers,paste("P33_",cluster.n,"clusters_all_markers.txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)

###################################################
####  cell type defination : ####
###################################################

DimPlot(s.object, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 7) +ggtitle("P33_named")+NoLegend()
#HC
VlnPlot( s.object, features = c("Pvalb","Epyc","Calb1","Cabp2","Tecta","Dcn","Ebf1","Spock1","Rgcc","Clic6","Efhd1",
                                "Hmgb2","Pcolce2","Tubb5","Tuba1b","Atoh1","Pcdh15","Pou4f3"), 
         pt.size = 0, ncol =5, assay = "RNA",slot = "counts",log = T)+NoLegend()


VlnPlot( s.object, features = c("Fgf8", "Atp2a3", "Rprm", "Kcnj13","Pcp4", "Slc26a5","Calca","Ocm", "Veph1", "Cacng2" , "Strip2", "Calb1","Cib2","Ighm","Ikzf2"), 
        pt.size = 0, ncol =5, assay = "RNA",slot = "counts",log = T)+NoLegend()
VlnPlot(s.object, features = c("Gata3","Atoh1","Fgfr3","Pou4f3","Chrna9","Ccer2", "Cabp2", "Calb2", "Pvalb","Pcdh15"), 
        pt.size = 0, ncol =5, assay = "RNA",slot = "counts",log = T)+NoLegend()
VlnPlot(s.object, features = c("Myo7a","Pou4f3","Gfi1","Pvalb","Otof","Calb2","Ptprq","Myo6","Dll1","Atoh1","Pax2","Dach1","Jag2","Gfi1","Slc17a8","Barhl1"), 
        pt.size = 0, ncol =5, assay = "RNA",slot = "counts",log = T)+NoLegend()

# -- candidate genes feature plot
#SC
VlnPlot(s.object, features = c("Fst","S100a6","Emid1","Npy","Sox2","S100a6","Pla2g7","Tuba1b","Spry2","Pcdh15"), 
        pt.size = 0, ncol =5, assay = "RNA",slot = "counts",log = T)+NoLegend()

VlnPlot(s.object, features = c("Jag1","Gli3","Notch3","Notch2","Maml2","Hes1","Hey1","Slc1a3","Egr3","Notch1","Sox9","Cdh1"), 
        pt.size = 0, ncol =5, assay = "RNA",slot = "counts",log = T)+NoLegend()
 
#DC
VlnPlot(s.object, features = c("Sox2","Rflnb","S100a1","Rbp7","Lfng","Hes5","Fgfr3"), 
        pt.size = 0, ncol =5, assay = "RNA",slot = "counts",log = T)+NoLegend()


P33_s.object.markers <- FindAllMarkers(s.object, min.pct = 0.25, logfc.threshold = 0.25)
write.table(P33_s.object.markers,"P33_cca.all_markers.xls",sep = "\t",quote = F,row.names = F)

# -- intersection plot
i=0  # select cluster id
cluster.markers <- s.object.markers[s.object.markers$cluster==i & s.object.markers$avg_log2FC>0,]
#top50_cluster0<-head(cluster0.markers, n = 50)
cluster.markers_0.05<-cluster.markers[cluster.markers$p_val_adj<=0.05 ,]

# hc : 9   (8)

# -- name cell types to clusters
cca.named<-readRDS("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P33\\P33_cca.renamed")
s.object<-cca.named
new.cluster.ids_2 <- c("Alas2+","HC","SC","SC","Cldn3+","SC","HC")
FeaturePlot(s.object, features = c("Atoh1"),
            pt.size =  0.05, ncol = 1)

names(new.cluster.ids_2) <- levels(s.object)
cca.named <- RenameIdents(s.object, new.cluster.ids_2)
DimPlot(cca.named, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) + NoLegend()
DimPlot(cca.named, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) + NoLegend()
saveRDS(cca.named, file = "P33_cca.renamed")

cca.named<-readRDS("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P33\\P33_cca.renamed")
##################P7_select_only_sc_hc########################
P33_select_SC_HC <- cca.named[,cca.named@meta.data$seurat_clusters %in% c(1,2,3,5,6)]
saveRDS(P33_select_SC_HC , file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select_SC_HC-seurat_foldchange\\P33_select_SC_HC")



##################P7_select_only_sc########################
P33_select_SC <- cca.named[,cca.named@meta.data$seurat_clusters %in% c(2,3,5)]
saveRDS(P33_select_SC , file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only SC\\foldchange\\P33_select_SC")

write.table(P33_select@assays$RNA@counts, "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only SC\\P33_select_only_SC_matrix.txt",sep="\t", quote=F, col.names=NA)



##################P7_select_only_hc########################

P33_select_hc <- cca.named[,cca.named@meta.data$seurat_clusters %in% c(1,6)]
saveRDS(P33_select_hc , file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only HC\\foldchange\\P33_select_hc")
write.table(P33_select_hc@assays$RNA@counts, "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only HC\\P33_select_only_HC_matrix.txt",sep="\t", quote=F, col.names=NA)

####################################################################################
#######################################################################################

#install.packages("clustree")
library(clustree)
clustree(s.object,prefix = "integrated_snn_res.")


install.packages("ggalluvial")
library(ggalluvial)
library(tidyverse)

ggplot(s.object@meta.data,aes(axis1 = integrated_snn_res.0.2, axis2 = integrated_snn_res.0.3,axis3 = integrated_snn_res.0.4,axis4 =integrated_snn_res.0.5,axis5 = integrated_snn_res.0.6,)) +
  scale_x_discrete(limits = c( "integrated_snn_res.0.2", "integrated_snn_res.0.3", "integrated_snn_res.0.4", "integrated_snn_res.0.5", "integrated_snn_res.0.6"))+
  geom_alluvium(aes(fill = integrated_snn_res.0.6), width= 1/12, show.legend = TRUE) +
  scale_colour_gradientn(colours=rainbow(17),guide = "colourbar") + 
  theme(legend.title = element_blank())+ 
  geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum))) + theme_bw()




