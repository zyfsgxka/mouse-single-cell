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
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P12")
P12_1.d <- Read10X(data.dir = "P12_1")
P12_2.d <- Read10X(data.dir = "P12_2")


# -- Create seurat object (non-normalized data)
## --- features: detected in > 3 cells
## --- cells: > 200 features are detected
P12_1<-CreateSeuratObject(count= P12_1.d, project = "P12_1", min.cells = 3, min.features = 200)
P12_2<-CreateSeuratObject(count= P12_2.d , project = "P12_2", min.cells = 3, min.features = 200)

# -- QC checking:
## --- calculate mitochondrial percent:
P12_1[["percent.mt"]] <- PercentageFeatureSet(P12_1, pattern = "^mt-")
P12_2[["percent.mt"]] <- PercentageFeatureSet(P12_2, pattern = "^mt-")

## --- label samples by development time:
P12_1[["time"]]  <-   "P12_1"
P12_2[["time"]]  <-   "P12_2"


## --- merge datasets: (NOTE: key step)
P12_merge <- merge(P12_1, y = c(P12_2), project = "P12_merge")
P12_merge # check cells

## --- Visualize QC metrics as a violin plot:
VlnPlot(P12_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(P12_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(P12_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# -- QC, clean datasets: (NOTE: key step, change threshold accordingly)
Rmin <- 200
Rmax <- 5000
mt.per <-5
P12_merge <- subset(P12_merge, subset = nFeature_RNA > Rmin & nFeature_RNA < Rmax & percent.mt < mt.per)
P12_merge # check cells


#################################################################
#### analysis for directly merged datasets: (for comparison) ####
#################################################################

s.object <- P12_merge  # define Seurat object

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
saveRDS(s.object,"P12_s.object")
## --- Examine and visualize PCA results a few different ways (optional)
print(s.object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(s.object, dims = 1:5, reduction = "pca")
DimHeatmap(s.object, dims = 1:15, cells = 500, balanced = TRUE) #cells = 500:seleceted 500 cells to run pca

## --- show PCA scores:  for determining the ??dimensionality??
npca <- 50  # define number of PCA scores
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
ndim <-30  # dims according to 'significant' PCs
## --- Computes the k.param nearest neighbors for a given dataset
s.object <- FindNeighbors(s.object, dims = 1:ndim)
## --- change resolution to get different clusters. (between 0.4-1.2 typically)
s.object <- FindClusters(s.object, resolution = 0.5) 
head(Idents(s.object), 5) # check cluster IDs of the first 5 cells

# -- Run non-linear dimensional reduction 
s.object <- RunUMAP(s.object, dims = 1:ndim)
p1 <- DimPlot(s.object, reduction = "umap",label = TRUE)
p1
p1_time <- DimPlot(s.object, reduction = "umap",label = TRUE,group.by = "time")
p1_time
p1+p1_time
saveRDS(s.object, file = "P12_merge.object")

###################################################
####  analysis for integrated datasets (CCA) : ####
###################################################

# -- split the dataset into a list of seurat objects
P12.cca <- SplitObject(P12_merge, split.by = "time")
#P1.cca <- SplitObject(s.object, split.by = "time")

# -- normalize and identify variable features for each dataset independently:
nG <- 2000
for (i in 1:length(P12.cca)){
  P12.cca[[i]]<- NormalizeData(P12.cca[[i]],verbose = FALSE)
  P12.cca[[i]]<- FindVariableFeatures(P12.cca[[i]],selection.method = "vst",
                                     nfeatures = nG, verbose = FALSE)
}

# -- Perform integration
reference.list <- P12.cca[c("P12_1","P12_2")]
CCA.anchors <- FindIntegrationAnchors(object.list=reference.list,dims = 1:30) 
## --- create an integrated data assay
P12.integrated <-IntegrateData(anchorset=CCA.anchors,dims = 1:30)

### -- analysis for integrated datasets -- ### 
# -- switch to integrated assay.
cca.object <- P12.integrated # define object
DefaultAssay(cca.object)<-"integrated" 

# -- Scaling the data
cca.object<-ScaleData(cca.object,verbose = FALSE)

# -- Perform linear dimensional reduction
cca.object<-RunPCA(cca.object,npcs = 50,verbose = FALSE)
print(cca.object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(cca.object, dims = 1:5, reduction = "pca")

npca <- 50  # number of PC scores to be calculated
## - identify ??significant?? PCs:
## opt-1:
ElbowPlot(cca.object,ndims = npca)
### opt-2: NOTE: The following process can take a LONG ~~ time
cca.object <- JackStraw(cca.object,dims = npca, num.replicate = 100)
cca.object <- ScoreJackStraw(cca.object, dims = 1:npca) 
JackStrawPlot(cca.object, dims = 1:npca)
#optional: output PCA done object 
#saveRDS(cca.object, file = "cca_PCAdone.object")
#cca.object <- readRDS("cca_PCAdone.object")

# -- Cluster the cells  
ndim <- 30  # dims according to 'significant' PCs
## --- Computes the k.param nearest neighbors for a given dataset
cca.object <- FindNeighbors(cca.object, dims = 1:ndim) 
## --- change resolution to get different clusters. (between 0.4-1.2 typically)

cca.object<- FindClusters(cca.object, resolution = 0.5)
head(cca.object)
# -- Run non-linear dimensional reduction
cca.object<-RunUMAP(cca.object,reduction = "pca",dims=1:ndim )
p3<-DimPlot(cca.object, reduction = "umap",group.by = "time")
p3
p4<-DimPlot(cca.object, reduction = "umap",label = TRUE,repel = TRUE)
p4
p3+p4
DimPlot(cca.object, reduction = "umap", split.by = "time") # visualize conditions side-by-side
#optional: output clustered object 
#saveRDS(cca.object, file = "P1_cca_clustered.object")
#cca.object <- readRDS("P1_cca_clustered.object")

# -- find markers for every cluster compared to all remaining cells (postive & negtive)
cca.object.markers <- FindAllMarkers(cca.object, min.pct = 0.25, logfc.threshold = 0.25)
#memory.limit()

#memory.limit(1000000)

# -- plotting the top 20 markers (or all markers if less than 20) for each cluster
top20 <- cca.object.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
DoHeatmap(cca.object, features = top20$gene)+scale_fill_gradientn(colors = c("skyblue", "white", "#FF3300"))

# -- output makers
cluster.n <- length(levels(cca.object.markers$cluster)) # get number of clusters
write.table(cca.object.markers,paste("cca.object_",cluster.n,"clusters_all_markers.txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)
saveRDS(cca.object, file = "P12.cca")


###################################################
####  cell type defination : ####
###################################################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P12")
cca.object<-readRDS('P12_cca.renamed')
DimPlot(cca.object, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 7) +ggtitle("P12_named")+NoLegend()

FeaturePlot(cca.object, features = c("Atoh1"),
            pt.size =  0.05, ncol = 1)


#HC
VlnPlot(cca.object, features = c("Fgf8", "Atp2a3", "Rprm", "Kcnj13","Pcp4", "Slc26a5","Calca","Ocm", "Veph1", "Cacng2" , "Strip2", "Calb1","Cib2","Ighm","Ikzf2"), 
        pt.size = 0, ncol =5, assay = "RNA",slot = "counts",log = T)+NoLegend()
VlnPlot(cca.object, features = c("Gata3","Atoh1","Fgfr3","Pou4f3","Chrna9","Ccer2", "Cabp2", "Calb2", "Pvalb"), 
        pt.size = 0, ncol =5, assay = "RNA",slot = "counts",log = T)+NoLegend()
# -- candidate genes feature plot
FeaturePlot(cca.object, features = c("Gata3","Atoh1","HA","Fgfr3","Pou4f3","Chrna9","Ccer2", "Acbd7", "Rprm", "Cd164l2"),
            pt.size =  0, ncol = 3,slot = "counts")+NoLegend()
#SC
VlnPlot(cca.object, features = c("Fst","S100a1","Emid1","Npy","Sox2"), 
        pt.size = 0, ncol =5, assay = "RNA",slot = "counts",log = T)+NoLegend()
#DC
VlnPlot(cca.object, features = c("Selenom","Gata3","Rbp7","Mansc4","Car13"), 
        pt.size = 0, ncol =5, assay = "RNA",slot = "counts",log = T)+NoLegend()


P12_cca.object.markers <- FindAllMarkers(cca.object, min.pct = 0.25, logfc.threshold = 0.25)
write.table(P12_cca.object.markers,"P12_cca.all_markers.xls",sep = "\t",quote = F,row.names = F)

# -- intersection plot
i=0  # select cluster id
cluster.markers <- cca.object.markers[cca.object.markers$cluster==i & cca.object.markers$avg_log2FC>0,]
#top50_cluster0<-head(cluster0.markers, n = 50)
cluster.markers_0.05<-cluster.markers[cluster.markers$p_val_adj<=0.05 ,]

# hc : 9   (8)

# -- name cell types to clusters
cca.object<-readRDS("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P12\\P12.cca")
new.cluster.ids_2 <- c("Notum+","SC","Slfn14+","Sdpr+","Itgam+","Opcml+","HC",
                       "HC","HC","SC","Bmp8a+","Mxd3+")
names(new.cluster.ids_2) <- levels(cca.object)
cca.named <- RenameIdents(cca.object, new.cluster.ids_2)

p1 <- DimPlot(cca.object, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) + NoLegend()
p2<-DimPlot(cca.named, reduction = "umap", label = TRUE,pt.size = 0.5,label.size =7) +ggtitle("P12")+ NoLegend()

p1 + p2
p2
cca.object<-readRDS("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P12\\P12_cca.renamed")
##################P12_select_only_sc########################
P12_select_SC_HC <- cca.object[,cca.object@meta.data$seurat_clusters %in% c(1,6,7,8,9)]
saveRDS(P12_select_SC_HC , file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select_SC_HC-seurat_foldchange\\P12_select_SC_HC")



##################P12_select_only_sc########################
P12_select_SC <- cca.object[,cca.object@meta.data$seurat_clusters %in% c(1,9)]
saveRDS(P12_select_SC , file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only SC\\foldchange\\P12_select_SC")
write.table(P12_select@assays$RNA@counts, "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only SC\\P12_select_only_SC_matrix.txt",sep="\t", quote=F, col.names=NA)


saveRDS(cca.named, file = "P12_cca.renamed")



##################P7_select_only_hc########################

P12_select_hc <- P12_cca.renamed[,P12_cca.renamed@meta.data$seurat_clusters %in% c(6,7,8)]

saveRDS(P12_select_hc , file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only HC\\foldchange\\P12_select_hc")
write.table(P12_select_hc@assays$RNA@counts, "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only HC\\P12_select_only_HC_matrix.txt",sep="\t", quote=F, col.names=NA)


####################################################################################
#######################################################################################

#install.packages("clustree")
library(clustree)
clustree(cca.object,prefix = "integrated_snn_res.")


install.packages("ggalluvial")
library(ggalluvial)
library(tidyverse)

ggplot(cca.object@meta.data,aes(axis1 = integrated_snn_res.0.2, axis2 = integrated_snn_res.0.3,axis3 = integrated_snn_res.0.4,axis4 =integrated_snn_res.0.5,axis5 = integrated_snn_res.0.6,)) +
  scale_x_discrete(limits = c( "integrated_snn_res.0.2", "integrated_snn_res.0.3", "integrated_snn_res.0.4", "integrated_snn_res.0.5", "integrated_snn_res.0.6"))+
  geom_alluvium(aes(fill = integrated_snn_res.0.6), width= 1/12, show.legend = TRUE) +
  scale_colour_gradientn(colours=rainbow(17),guide = "colourbar") + 
  theme(legend.title = element_blank())+ 
  geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum))) + theme_bw()




