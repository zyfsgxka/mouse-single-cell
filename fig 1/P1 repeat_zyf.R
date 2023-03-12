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
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P1_repeat")
P1_1.d <- Read10X(data.dir = "p1_1")
P1_2.d <- Read10X(data.dir = "p1_2")
P1_3.d <- Read10X(data.dir = "p1_3")
P1_4.d <- Read10X(data.dir = "p1_4")

# -- Create seurat object (non-normalized data)
## --- features: detected in > 3 cells
## --- cells: > 200 features are detected
P1_1<-CreateSeuratObject(count= P1_1.d, project = "P1_1", min.cells = 3, min.features = 200)
P1_2<-CreateSeuratObject(count= P1_2.d , project = "P1_2", min.cells = 3, min.features = 200)
P1_3<-CreateSeuratObject(count= P1_3.d , project = "P1_3", min.cells = 3, min.features = 200)
P1_4<-CreateSeuratObject(count= P1_4.d , project = "P1_4", min.cells = 3, min.features = 200)

# -- QC checking:
## --- calculate mitochondrial percent:
P1_1[["percent.mt"]] <- PercentageFeatureSet(P1_1, pattern = "^mt-")
P1_2[["percent.mt"]] <- PercentageFeatureSet(P1_2, pattern = "^mt-")
P1_3[["percent.mt"]] <- PercentageFeatureSet(P1_3, pattern = "^mt-")
P1_4[["percent.mt"]] <- PercentageFeatureSet(P1_4, pattern = "^mt-")
## --- label samples by development time:
P1_1[["time"]]  <-   "P1_1"
P1_2[["time"]]  <-   "P1_2"
P1_3[["time"]]  <-   "P1_3"
P1_4[["time"]]  <-   "P1_4"

## --- merge datasets: (NOTE: key step)
P1_merge <- merge(P1_1, y = c(P1_2,P1_3,P1_4), project = "P1_merge")
P1_merge # check cells

## --- Visualize QC metrics as a violin plot:
VlnPlot(P1_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(P1_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(P1_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# -- QC, clean datasets: (NOTE: key step, change threshold accordingly)
Rmin <- 200
Rmax <- 4500
mt.per <-5
P1_merge <- subset(P1_merge, subset = nFeature_RNA > Rmin & nFeature_RNA < Rmax & percent.mt < mt.per)
P1_merge # check cells


#################################################################
#### analysis for directly merged datasets: (for comparison) ####
#################################################################

s.object <- P1_merge  # define Seurat object

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
s.object <- FindClusters(s.object, resolution = 0.2) 
head(Idents(s.object), 5) # check cluster IDs of the first 5 cells

# -- Run non-linear dimensional reduction 
s.object <- RunUMAP(s.object, dims = 1:ndim)
p1 <- DimPlot(s.object, reduction = "umap",label = TRUE)
p1
p1_time <- DimPlot(s.object, reduction = "umap",label = TRUE,group.by = "time")
p1_time
p1+p1_time
saveRDS(s.object, file = "P1_s.object")

###################################################
####  analysis for integrated datasets (CCA) : ####
###################################################

# -- split the dataset into a list of seurat objects
P1.cca <- SplitObject(P1_merge, split.by = "time")
#P1.cca <- SplitObject(s.object, split.by = "time")

# -- normalize and identify variable features for each dataset independently:
nG <- 2000
for (i in 1:length(P1.cca)){
  P1.cca[[i]]<- NormalizeData(P1.cca[[i]],verbose = FALSE)
  P1.cca[[i]]<- FindVariableFeatures(P1.cca[[i]],selection.method = "vst",
                                      nfeatures = nG, verbose = FALSE)
}

# -- Perform integration
reference.list <- P1.cca[c("P1_1","P1_2","P1_3","P1_4")]
CCA.anchors <- FindIntegrationAnchors(object.list=reference.list,dims = 1:30) 
## --- create an integrated data assay
P1.integrated <-IntegrateData(anchorset=CCA.anchors,dims = 1:30)

### -- analysis for integrated datasets -- ### 
# -- switch to integrated assay.
cca.object <- P1.integrated # define object
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
ndim <- 40  # dims according to 'significant' PCs
## --- Computes the k.param nearest neighbors for a given dataset
cca.object <- FindNeighbors(cca.object, dims = 1:ndim) 
## --- change resolution to get different clusters. (between 0.4-1.2 typically)

cca.object<- FindClusters(cca.object, resolution = 0.5)


head(cca.object)

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



memory.limit()

memory.limit(1000000)

# -- plotting the top 20 markers (or all markers if less than 20) for each cluster
top20 <- cca.object.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
DoHeatmap(cca.object, features = top20$gene)+scale_fill_gradientn(colors = c("skyblue", "white", "#FF3300"))

# -- output makers
cluster.n <- length(levels(cca.object.markers$cluster)) # get number of clusters
write.table(cca.object.markers,paste("cca.object_",cluster.n,"clusters_all_markers.txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)

###################################################
####  cell type defination : ####
###################################################

# -- input reference genes from paper
ref_marker.d <- read.csv("P1_top_gene_paper.txt",sep="\t",header = T)
cell_type <- names(ref_marker.d)
names(ref_marker.d)  # ref cell types

# -- check marker genes expression (violin plot by clusters)
for(i in 1:length(cell_type)){
  cell_type[i]  # check cell type
  genes <- ref_marker.d[1:20,cell_type[i]]
  pdf(paste(cell_type[i],".pdf",sep=""),width = 20,height =14 )
  p <- VlnPlot(cca.object,assay="RNA",features = genes,pt.size = 0)
  print(p)
  dev.off()
}
saveRDS(cca.object, file = "P1.cca")

# -- candidate genes feature plot
FeaturePlot(cca.object, features = c("Ccer2","Selenom","Pcp4","Lhfpl5","Pou4f3","Cib2","Evl","Lmo1","Pvalb","Grxcr1","Atoh1","Grp","Hes6","Psph","Gng8","Acbd7","Myo6","Ccl21a","Gm2694","Calm1"),
            pt.size =  0, ncol = 5)

# -- intersection plot
i=0  # select cluster id
cluster.markers <- cca.object.markers[cca.object.markers$cluster==i & cca.object.markers$avg_log2FC>0,]
#top50_cluster0<-head(cluster0.markers, n = 50)
cluster.markers_0.05<-cluster.markers[cluster.markers$p_val_adj<=0.05 ,]

library(UpSetR)
input <- list(maker=cluster.markers$gene,
              gene_OHC=ref_marker.d$gene_OHC,
              gene_KO=ref_marker.d$gene_KO,
              gene_Hensen=ref_marker.d$gene_Hensen,
              gene_IHC=ref_marker.d$gene_IHC,
              gene_Interdental=ref_marker.d$gene_Interdental,
              gene_IPC=ref_marker.d$gene_IPC,
              gene_IS=ref_marker.d$gene_IS,
              gene_Lateral_Prosensory=ref_marker.d$gene_Lateral_Prosensory,
              gene_Lateral_KO=ref_marker.d$gene_Lateral_KO,
              gene_LER_Bmp4=ref_marker.d$gene_LER_Bmp4,
              gene_LER_Fst=ref_marker.d$gene_LER_Fst,
              gene_Medial_Prosensory=ref_marker.d$gene_Medial_Prosensory,
              gene_Oc90_Otoa=ref_marker.d$gene_Oc90_Otoa,
              gene_Oc90_Sparcl1=ref_marker.d$gene_Oc90_Sparcl1,
              gene_iOHC=ref_marker.d$gene_iOHC,
              gene_OS_Claudius=ref_marker.d$gene_OS_Claudius,
              gene_Ube2c=ref_marker.d$gene_Ube2c)

upset(fromList(input),nsets = 20,)

# -- name cell types to clusters
new.cluster.ids_2 <- c("M.KO","L.KO","L.KO","unknown","unknown","ISC/IdC","unknown","unknown","unknown","HC",
                       "DC/OPC","unknown","CC/OSC/HeC","IPhc","unknown","unknown","OC90","unknown","unknown","HC","IPC","unknown","unknown","unknown")
names(new.cluster.ids_2) <- levels(cca.object)
cca.named <- RenameIdents(cca.object, new.cluster.ids_2)
DimPlot(cca.named, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) + NoLegend()

saveRDS(cca.named, file = "P1.cca.named")


######rename

setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2\\P1")
P1_cca<-readRDS("P1.cca")
DimPlot(P1_cca, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 7) +ggtitle("P1_named")+NoLegend()

P1_cca.markers <- FindAllMarkers(P1_cca, min.pct = 0.25, logfc.threshold = 0.25)

write.table(P1_cca.markers,"P1_cca.all_markers.xls",sep = "\t",quote = F,row.names = F)

new.cluster.ids_3 <- c("M.KO","L.KO","L.KO","Sptssa+","Mgp+","ISC/IdC","Otor+",
                       "Otor","Ccn3+","HC","DC/OPC","Arpc1b+","CC/OSC/HeC",
                       "IPhc","Adgre1+","Hmgb2+",
                       "OC90","S100a1+","Gm14635+","HC","IPC","Trpm1+","Cpne4+",
                       "Esrrb+")


names(new.cluster.ids_3) <- levels(P1_cca)
cca.renamed <- RenameIdents(P1_cca, new.cluster.ids_3)
p5<-DimPlot(cca.renamed, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 5) +ggtitle("P1_named")+NoLegend()
p4+p5
p5

saveRDS(cca.renamed,file ="P1_cca.renamed")



##################P1_select_sc_hc########################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P1_repeat")
P1_cca.renamed<-readRDS('P1_cca.renamed')
P1_select_SC_HC <- P1_cca.renamed[,P1_cca.renamed@meta.data$seurat_clusters %in% c(9,10,12,13,19,20)]
saveRDS(P1_select_SC_HC , file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select_SC_HC-seurat_foldchange\\P1_select_SC_HC")

write.table(P1_select@assays$RNA@counts, "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2\\P1_select\\P1_select_matrix.txt",sep="\t", quote=F, col.names=NA)

##################P1_select_sc########################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P1_repeat")
P1_cca.renamed<-readRDS('P1_cca.renamed')
P1_select_SC <- P1_cca.renamed[,P1_cca.renamed@meta.data$seurat_clusters %in% c(10,12,13,20)]
saveRDS(P1_select_SC , file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only SC\\foldchange\\P1_select_SC")
write.table(P1_select@assays$RNA@counts, "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only SC\\P1_select_only_SC_matrix.txt",sep="\t", quote=F, col.names=NA)

FeaturePlot(P1_cca.renamed, features = c("Atoh1"),assay = "RNA",slot = "counts",
            pt.size =  0.05, ncol = 1)


##################P1_select_hc########################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P1_repeat")
P1_cca.renamed<-readRDS('P1_cca.renamed')
P1_select_hc <- P1_cca.renamed[,P1_cca.renamed@meta.data$seurat_clusters %in% c(9,19)]
saveRDS(P1_select_hc , file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only HC\\foldchange\\P1_select_hc")
write.table(P1_select_hc@assays$RNA@counts, "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only HC\\P1_select_only_HC_matrix.txt",sep="\t", quote=F, col.names=NA)
