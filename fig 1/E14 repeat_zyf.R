library(Seurat)
library(dplyr)
library(magrittr)
library(hdf5r)
library(ggplot2)
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\E14 repeat")
#E14_1 <- Read10X(data.dir = "E14_1")


E14_2<-readRDS('E14_2')

data("E14_2")
d<-FoldChange(E14_2, ident.1 = 1)

d


E14_2 <- Read10X(data.dir = "E14_2 primary files")
##create seurat object
#数据集中测到的少于200个基因的细胞（min.features = 200）和少于3个细胞覆盖的基因（min.cells = 3）被过滤掉
#E14_1<-CreateSeuratObject(count= E14_1, project = "E14_1", min.cells = 3, min.features = 200)
E14_2<-CreateSeuratObject(count= E14_2 , project = "E14_2", min.cells = 3, min.features = 200)



#E14_1[["percent.mt"]] <- PercentageFeatureSet(E14_1, pattern = "^mt-")
E14_2[["percent.mt"]] <- PercentageFeatureSet(E14_2, pattern = "^mt-")


#E14_1[["time"]]  <-   "E14_1"
E14_2[["time"]]  <-   "E14_2"



#E14 and E16 merge assay
#E14_1_2_repeat_merge <- merge(E14_2, y = c(E14_1), project = "E14_1_2_repeat_merge")
VlnPlot(E14_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(E14_2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E14_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
E14_2<- subset(E14_2, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 5)
E14_2



E14_2 <- NormalizeData(E14_2, normalization.method = "LogNormalize", scale.factor = 10000)
E14_2 <- FindVariableFeatures(E14_2, selection.method = "vst", nfeatures = 2000)#nfeatures = 2000:	Number of features to select as top variable features
top20<- head(VariableFeatures(E14_2), 20)
top20

plot1 <- VariableFeaturePlot(E14_2)
plot2 <- LabelPoints(plot = plot1, points = top20)
plot1 + plot2
plot2
#plot2+theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),
#            axis.title.x= element_text(size=15),axis.title.y= element_text(size=15)


all.genes <- rownames(E14_2)
E14_2<- ScaleData(E14_2, features = all.genes)

E14_2 <- RunPCA(E14_2, features = VariableFeatures(object = E14_2))
#print(E14_2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(E14_2, dims = 1:2, reduction = "pca")
DimPlot(E14_2, reduction = "pca")
#dim:	Dimensions to plot,indicated to see which pca; cells = 500:seleceted 500 cells to run pca
DimHeatmap(E14_2, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(E14_2, dims = 1:15, cells = 500, balanced = TRUE)

#num.replicate:Number of replicate samplings to perform
#Randomly permutes a subset of data, and calculates projected PCA scores for these 'random' genes. 
#Then compares the PCA scores for the 'random' genes with the observed PCA scores to determine statistical signifance. 
#End result is a p-value for each gene's association with each principal component.
E14_2 <- JackStraw(E14_2, num.replicate = 100)
E14_2 <- ScoreJackStraw(E14_2, dims = 1:20)
JackStrawPlot(E14_2, dims = 1:20)
ElbowPlot(E14_2)

#Computes the k.param nearest neighbors for a given dataset.
E14_2 <- FindNeighbors(E14_2, dims = 1:20)
E14_2 <- FindClusters(E14_2, resolution = 0.5)
head(Idents(E14_2), 5)

E14_2 <- RunUMAP(E14_2, dims = 1:20)
p1_time <- DimPlot(E14_2, reduction = "umap",label = TRUE,group.by = "time")
p1 <- DimPlot(E14_2, reduction = "umap",label = TRUE)
p1_time
p1


  
E14_2.markers <- FindAllMarkers(E14_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
E14_2.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
top20 <- E14_2.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(E14_2, features = top20$gene) + NoLegend()


saveRDS(E14_2, file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\重复repeat\\E14 repeat\\E14_2")


P=0.05

cluster0.markers <- FindMarkers(E14_2, ident.1 = 0, min.pct = 0.25)

#top50_cluster0<-head(cluster0.markers, n = 50)
cluster0.markers_0.05<-cluster0.markers[cluster0.markers$p_val_adj<=P,]
write.table(cluster0.markers_0.05,"cluster0.markers",sep = "\t",quote = F,row.names = T,col.names = T)
cluster0.markers<- read.csv("cluster0.markers",sep = "\t")
cluster0.markers$gene<- row.names(cluster0.markers)

cluster1.markers <- FindMarkers(E14_2, ident.1 = 1, min.pct = 0.25)
#top50_cluster1<-head(cluster1.markers, n = 50)
cluster1.markers_0.05<-cluster1.markers[cluster1.markers$p_val_adj<=P,]
write.table(cluster1.markers_0.05,"cluster1.markers",sep = "\t",quote = F,row.names = T,col.names = T)
cluster1.markers<- read.csv("cluster1.markers",sep = "\t")
cluster1.markers$gene<- row.names(cluster1.markers)

cluster2.markers <- FindMarkers(E14_2, ident.1 = 2, min.pct = 0.25)
#top50_cluster2<-head(cluster2.markers, n = 50)
cluster2.markers_0.05<-cluster2.markers[cluster2.markers$p_val_adj<=P,]
write.table(cluster2.markers_0.05,"cluster2.markers",sep = "\t",quote = F,row.names = T,col.names = T)
cluster2.markers<- read.csv("cluster2.markers",sep = "\t")
cluster2.markers$gene<- row.names(cluster2.markers)

cluster3.markers <- FindMarkers(E14_2, ident.1 = 3, min.pct = 0.25)
#top50_cluster3<-head(cluster3.markers, n = 50)
cluster3.markers_0.05<-cluster3.markers[cluster3.markers$p_val_adj<=P,]
write.table(cluster3.markers_0.05,"cluster3.markers",sep = "\t",quote = F,row.names = T,col.names = T)
cluster3.markers<- read.csv("cluster3.markers",sep = "\t")
cluster3.markers$gene<- row.names(cluster3.markers)

cluster4.markers <- FindMarkers(E14_2, ident.1 = 4, min.pct = 0.25)
#top50_cluster4<-head(cluster4.markers, n = 50)
cluster4.markers_0.05<-cluster4.markers[cluster4.markers$p_val_adj<=P,]
write.table(cluster4.markers_0.05,"cluster4.markers",sep = "\t",quote = F,row.names = T,col.names = T)
cluster4.markers<- read.csv("cluster4.markers",sep = "\t")
cluster4.markers$gene<- row.names(cluster4.markers)

cluster5.markers <- FindMarkers(E14_2, ident.1 = 5, min.pct = 0.25)
#top50_cluster5<-head(cluster5.markers, n = 50)
cluster5.markers_0.05<-cluster5.markers[cluster5.markers$p_val_adj<=P,]
write.table(cluster5.markers_0.05,"cluster5.markers",sep = "\t",quote = F,row.names = T,col.names = T)
cluster5.markers<- read.csv("cluster5.markers",sep = "\t")
cluster5.markers$gene<- row.names(cluster5.markers)

cluster6.markers <- FindMarkers(E14_2, ident.1 = 6, min.pct = 0.25)
#top50_cluster6<-head(cluster6.markers, n = 50)
cluster6.markers_0.05<-cluster6.markers[cluster6.markers$p_val_adj<=P,]
write.table(cluster6.markers_0.05,"cluster6.markers",sep = "\t",quote = F,row.names = T,col.names = T)
cluster6.markers<- read.csv("cluster6.markers",sep = "\t")
cluster6.markers$gene<- row.names(cluster6.markers)


cluster7.markers <- FindMarkers(E14_2, ident.1 = 7, min.pct = 0.25)
#top50_cluster7<-head(cluster7.markers, n = 50)
cluster7.markers_0.05<-cluster7.markers[cluster7.markers$p_val_adj<=P,]
write.table(cluster7.markers_0.05,"cluster7.markers",sep = "\t",quote = F,row.names = T,col.names = T)
cluster7.markers<- read.csv("cluster7.markers",sep = "\t")
cluster7.markers$gene<- row.names(cluster7.markers)

cluster8.markers <- FindMarkers(E14_2, ident.1 = 8, min.pct = 0.25)
#top50_cluster8<-head(cluster8.markers, n = 50)
cluster8.markers_0.05<-cluster8.markers[cluster8.markers$p_val_adj<=P,]
write.table(cluster8.markers_0.05,"cluster8.markers ",sep = "\t",quote = F,row.names = T,col.names = T)
cluster8.markers<- read.csv("cluster8.markers",sep = "\t")
cluster8.markers$gene<- row.names(cluster8.markers)

cluster9.markers <- FindMarkers(E14_2, ident.1 = 9, min.pct = 0.25)
#top50_cluster9<-head(cluster9.markers, n = 50)
cluster9.markers_0.05<-cluster9.markers[cluster9.markers$p_val_adj<=P,]
write.table(cluster9.markers_0.05,"cluster9.markers",sep = "\t",quote = F,row.names = T,col.names = T)
cluster9.markers<- read.csv("cluster9.markers",sep = "\t")
cluster9.markers$gene<- row.names(cluster9.markers)


cluster10.markers <- FindMarkers(E14_2, ident.1 = 10, min.pct = 0.25)
#top50_cluster10<-head(cluster10.markers, n = 50)
cluster10.markers_0.05<-cluster10.markers[cluster10.markers$p_val_adj<=P,]
write.table(cluster10.markers_0.05,"cluster10.markers ",sep = "\t",quote = F,row.names = T,col.names = T)
cluster10.markers<- read.csv("cluster10.markers",sep = "\t")
cluster10.markers$gene<- row.names(cluster10.markers)

cluster11.markers <- FindMarkers(E14_2, ident.1 = 11, min.pct = 0.25)
#top50_cluster11<-head(cluster11.markers, n = 50)
cluster11.markers_0.05<-cluster11.markers[cluster11.markers$p_val_adj<=P,]
write.table(cluster11.markers_0.05,"cluster11.markers ",sep = "\t",quote = F,row.names = T,col.names = T)
cluster11.markers<- read.csv("cluster11.markers",sep = "\t")
cluster11.markers$gene<- row.names(cluster11.markers)

cluster12.markers <- FindMarkers(E14_2, ident.1 = 12, min.pct = 0.25)
#top50_cluster12<-head(cluster12.markers, n = 50)
cluster12.markers_0.05<-cluster12.markers[cluster12.markers$p_val_adj<=P,]
write.table(cluster12.markers_0.05,"cluster12.markers",sep = "\t",quote = F,row.names = T,col.names = T)
cluster12.markers<- read.csv("cluster12.markers",sep = "\t")
cluster12.markers$gene<- row.names(cluster12.markers)














E14_2_top_gene_paper<- read.table("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\重复repeat\\E14 repeat\\E14_2_top gene_paper.txt",sep = "\t",header = T)
E14_2_top_gene_paper$gene_KO <- row.names(E14_2_top_gene_paper_gene_KO)
E14_2_top_gene_paper$gene_KO/Hmgn2 <- row.names(E14_2_top_gene_paper_gene_KO/Hmgn2)
E14_2_top_gene_paper$gene_IHC <- row.names(E14_2_top_gene_paper_gene_IHC)
E14_2_top_gene_paper$gene_Interdental <- row.names(E14_2_top_gene_paper_gene_Interdental)
E14_2_top_gene_paper$gene_Lateral_Prosensory <- row.names(E14_2_top_gene_paper_gene_Lateral_Prosensory)
E14_2_top_gene_paper$gene_LER/Bmp4_1 <- row.names(E14_2_top_gene_paper_gene_LER/Bmp4_1)
E14_2_top_gene_paper$gene_LER/Fst <- row.names(E14_2_top_gene_paper_gene_LER/Fst)
E14_2_top_gene_paper$gene_Medial_Prosensory <- row.names(E14_2_top_gene_paper_gene_Medial_Prosensory)
E14_2_top_gene_paper$gene_Oc90_3 <- row.names(E14_2_top_gene_paper_gene_Oc90_3)
E14_2_top_gene_paper$gene_Oc90_4 <- row.names(E14_2_top_gene_paper_gene_Oc90_4)
E14_2_top_gene_paper$gene_Oc90/Otoa <- row.names(E14_2_top_gene_paper_gene_Oc90/Otoa)
E14_2_top_gene_paper$gene_Oc90/Sparcl1 <- row.names(E14_2_top_gene_paper_gene_Oc90/Sparcl1)
E14_2_top_gene_paper$gene_OHC <- row.names(E14_2_top_gene_paper_gene_OHC)

top50_cluster0<- read.csv("top50_cluster0",sep = "\t")
top50_cluster0$gene<- row.names(top50_cluster0)

top50_cluster1<- read.csv("top50_cluster1",sep = "\t")
top50_cluster1$gene<- row.names(top50_cluster1)

top50_cluster2<- read.csv("top50_cluster2",sep = "\t")
top50_cluster2$gene<- row.names(top50_cluster2)

top50_cluster3<- read.csv("top50_cluster3",sep = "\t")
top50_cluster3$gene<- row.names(top50_cluster3)

top50_cluster4<- read.csv("top50_cluster4",sep = "\t")
top50_cluster4$gene<- row.names(top50_cluster4)

top50_cluster5<- read.csv("top50_cluster5",sep = "\t")
top50_cluster5$gene<- row.names(top50_cluster5)

top50_cluster6<- read.csv("top50_cluster6",sep = "\t")
top50_cluster6$gene<- row.names(top50_cluster6)

top50_cluster7<- read.csv("top50_cluster7",sep = "\t")
top50_cluster7$gene<- row.names(top50_cluster7)

top50_cluster8<- read.csv("top50_cluster8",sep = "\t")
top50_cluster8$gene<- row.names(top50_cluster8)

top50_cluster9<- read.csv("top50_cluster9",sep = "\t")
top50_cluster9$gene<- row.names(top50_cluster9)

top50_cluster10<- read.csv("top50_cluster10",sep = "\t")
top50_cluster10$gene<- row.names(top50_cluster10)

top50_cluster11<- read.csv("top50_cluster11",sep = "\t")
top50_cluster11$gene<- row.names(top50_cluster11)

top50_cluster12<- read.csv("top50_cluster12",sep = "\t")
top50_cluster12$gene<- row.names(top50_cluster12)




#BiocManager::install("UpSetR")
library(UpSetR)
input <- list(cluster12=cluster12.markers$gene,
                  gene_KO=E14_2_top_gene_paper$gene_KO,
                  gene_KO.Hmgn2=E14_2_top_gene_paper$gene_KO.Hmgn2,
                  gene_IHC=E14_2_top_gene_paper$gene_IHC,
                  gene_Interdental=E14_2_top_gene_paper$gene_Interdental,
                  gene_Lateral.Prosensory=E14_2_top_gene_paper$gene_Lateral.Prosensory,
                  gene_LER.Bmp4_1=E14_2_top_gene_paper$gene_LER.Bmp4_1,
                  gene_LER.Fst=E14_2_top_gene_paper$gene_LER.Fst,
                  gene_Medial_Prosensory=E14_2_top_gene_paper$gene_Medial_Prosensory,
                  gene_Oc90_3=E14_2_top_gene_paper$gene_Oc90_3,
                  gene_Oc90_4=E14_2_top_gene_paper$gene_Oc90_4,
                  gene_Oc90.Otoa=E14_2_top_gene_paper$gene_Oc90.Otoa,
                  gene_Oc90.Sparcl1=E14_2_top_gene_paper$gene_Oc90.Sparcl1,
                  gene_OHC=E14_2_top_gene_paper$gene_OHC)


upset(fromList(input),nsets = 14,keep.order=TRUE)

BiocManager::install


#KO（1?,2,7）
VlnPlot(E14_2, features = c("Epyc","Calb1","Pclaf","Tecta","Dcn","Ebf1","1500015O10Rik","Spock1","Rgcc","Clic6","Efhd1","Hmgb2","Pcolce2","Tubb5","Tuba1b","Itga8","Snx3","Bcl2","H2afx","2410124H12Rik"), 
        pt.size = 0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Epyc","Calb1","Pclaf","Tecta","Dcn","Ebf1","1500015O10Rik","Spock1","Rgcc","Clic6","Efhd1","Hmgb2","Pcolce2","Tubb5","Tuba1b","Itga8","Snx3","Bcl2","H2afx","2410124H12Rik"),
            pt.size =  0, ncol = 5)

#KO/Hmgn2(2,6,10?)
VlnPlot(E14_2, features = c("Cenpa","Hist1h2bc","Pttg1","Ccnb2","Stmn1","Cdc20","H2afz","Hmgn2","2610318N02Rik","Hmgb2","Birc5","Hmgb1","Hist1h1c","Cdca8","H2afv","Cdca3","Knstrn","Acsl5","Hmgb3","Selenoh"), 
        pt.size = 0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Cenpa","Hist1h2bc","Pttg1","Ccnb2","Stmn1","Cdc20","H2afz","Hmgn2","2610318N02Rik","Hmgb2","Birc5","Hmgb1","Hist1h1c","Cdca8","H2afv","Cdca3","Knstrn","Acsl5","Hmgb3","Selenoh"),
            pt.size =  0.05, ncol = 5)

#IHC(11)
VlnPlot(E14_2, features = c("Ccer2","Pcp4","Selenom","Pou4f3","S100a1","Atoh1","Fgf8","Cib2","Hes6","Gm2694","Grxcr1","Lrrc73","Evl","Lmo1","Tmem255b","Calm1","Rprm","Rasd2","Dlk2","Psph"), 
        pt.size = 0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Ccer2","Pcp4","Selenom","Pou4f3","S100a1","Atoh1","Fgf8","Cib2","Hes6","Gm2694","Grxcr1","Lrrc73","Evl","Lmo1","Tmem255b","Calm1","Rprm","Rasd2","Dlk2","Psph"),
            pt.size =  0.05, ncol = 5)

#Interdental(1,2,7,8) similar with ko？
VlnPlot(E14_2, features = c("Smoc2","Foxq1","1500015O10Rik","Cdkn1c","Gsn","Tpm2","Ptn","Mia","Slc12a2","Otoa","Pla2r1","Ckb","Crip2","H19","Fbp2","Mfap4","Aldh1a3","Rgcc","Mettl9","Galm"), 
        pt.size = 0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Smoc2","Foxq1","1500015O10Rik","Cdkn1c","Gsn","Tpm2","Ptn","Mia","Slc12a2","Otoa","Pla2r1","Ckb","Crip2","H19","Fbp2","Mfap4","Aldh1a3","Rgcc","Mettl9","Galm"),
            pt.size =  0.05, ncol = 5)


#Lateral Prosensory（3.4.11）
VlnPlot(E14_2, features = c("Socs2","Fgfr3","Rprm","Camta1","Lockd","Mdm1","Ntf3","Hpcal1","Phlda1","Uchl1","Hes6","Hey1","S100a1","Mgst3","Isl1","Fkbp1a","Rtn1","Skp1a","Rcn1","Lor"), 
        pt.size =0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Socs2","Fgfr3","Rprm","Camta1","Lockd","Mdm1","Ntf3","Hpcal1","Phlda1","Uchl1","Hes6","Hey1","S100a1","Mgst3","Isl1","Fkbp1a","Rtn1","Skp1a","Rcn1","Lor"),
            pt.size =  0.05, ncol = 5)


#LER/Bmp4_1(4)
VlnPlot(E14_2, features = c("Bmp4","Fbln2","Gata2","Clu","Hs3st1","Plet1","Fst","C1qtnf12","Gata3","Lmo3","Itih5","Fbxo2","App","Tac1","Rtn1","Grb10","Sparc","Sfrp1","Serpine2","Id3"), 
        pt.size =0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Bmp4","Fbln2","Gata2","Clu","Hs3st1","Plet1","Fst","C1qtnf12","Gata3","Lmo3","Itih5","Fbxo2","App","Tac1","Rtn1","Grb10","Sparc","Sfrp1","Serpine2","Id3"),
            pt.size =  0.05, ncol = 5)

#LER/Fst（4,11?）
VlnPlot(E14_2, features = c("Tac1","Fst","Clu","Hmga2","Rprm","Isl1","Gata3","Hs3st1","Nr2f2","Igfbp5","App","Msx1","Plet1","Scara3","Emx2","Itih5","Angpt1","C1qtnf12","Nptx2","Vim"), 
        pt.size = 0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Tac1","Fst","Clu","Hmga2","Rprm","Isl1","Gata3","Hs3st1","Nr2f2","Igfbp5","App","Msx1","Plet1","Scara3","Emx2","Itih5","Angpt1","C1qtnf12","Nptx2","Vim"),
            pt.size =  0.05, ncol = 5)


#Medial_Prosensory（3,4,11）
VlnPlot(E14_2, features = c("Tectb","Socs2","S100a1","Anxa5","Isl1","Igfbp3","Lor","Fgf20","Ntf3","Scg5","Pcp4","Pcp4l1","Crym","Ntm","Hmga2","Phlda1","Moxd1","Cyp26b1","Skp1a","Sox2"), 
        pt.size = 0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Tectb","Socs2","S100a1","Anxa5","Isl1","Igfbp3","Lor","Fgf20","Ntf3","Scg5","Pcp4","Pcp4l1","Crym","Ntm","Hmga2","Phlda1","Moxd1","Cyp26b1","Skp1a","Sox2"),
            pt.size =  0.05, ncol = 5)

#Oc90_3（5，8)
VlnPlot(E14_2, features = c("Oc90","Ecel1","Cndp2","Kctd8","Fibin","Sptssa","Meg3","Dapl1","Krt19","Cyp26a1","Tbx1","Arl4c","Gpd2","BC006965","Cnmd","Ppp2r2b","Otx2","Bricd5","Entpd3","Wfdc2"), 
        pt.size = 0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Oc90","Ecel1","Cndp2","Kctd8","Fibin","Sptssa","Meg3","Dapl1","Krt19","Cyp26a1","Tbx1","Arl4c","Gpd2","BC006965","Cnmd","Ppp2r2b","Otx2","Bricd5","Entpd3","Wfdc2"),
            pt.size =  0.05, ncol = 5)

#Oc90_4(5,8)
VlnPlot(E14_2, features = c("Ecel1","Oc90","Igfbp5","Anxa2","Frem2","Anxa5","Id4","Arl4c","Gata2","Kctd8","Cebpb","Tbx1","Irx1","Plscr2","Sbspon","Maff","Lrpap1","Enpep","Fbn2","Dlx2"), 
        pt.size = 0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Ecel1","Oc90","Igfbp5","Anxa2","Frem2","Anxa5","Id4","Arl4c","Gata2","Kctd8","Cebpb","Tbx1","Irx1","Plscr2","Sbspon","Maff","Lrpap1","Enpep","Fbn2","Dlx2"),
            pt.size =  0.05, ncol = 5)


#Oc90/Otoa(8)
VlnPlot(E14_2, features = c("Ttr","Slc6a15","Oc90","Dapl1","Ptgds","Postn","Cndp2","Meg3","F3","Cacna2d1","Agr2","Krt19","Calml4","Sulf1","B4galt6","Cthrc1","Vmo1","Mme","Lurap1l","Ptn"), 
        pt.size = 0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Ttr","Slc6a15","Oc90","Dapl1","Ptgds","Postn","Cndp2","Meg3","F3","Cacna2d1","Agr2","Krt19","Calml4","Sulf1","B4galt6","Cthrc1","Vmo1","Mme","Lurap1l","Ptn"),
            pt.size =  0.05, ncol = 5)

#Oc90/Sparcl1(5)
VlnPlot(E14_2, features = c("Sptssa","Oc90","Ecel1","Anxa1","Lrpap1","Enpep","Meg3","Ccdc3","Gsc","Fam20a","Id4","Cldn4","Lrp2","Igf1","Vwc2","Kctd8","Aqp5","Sbspon","Washc3","Podxl2"), 
        pt.size = 0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Sptssa","Oc90","Ecel1","Anxa1","Lrpap1","Enpep","Meg3","Ccdc3","Gsc","Fam20a","Id4","Cldn4","Lrp2","Igf1","Vwc2","Kctd8","Aqp5","Sbspon","Washc3","Podxl2"),
            pt.size =  0.05, ncol = 5)


#OHC(11)
VlnPlot(E14_2, features = c("Ccer2","Pcp4","Pou4f3","Hes6","Rprm","Psph","Atoh1","Gadd45g","Calm1","Rasd2","Tesc","Fgfr3","Hpcal1","Grp","Tmem255b","Dlk2","Cryab","Kif19a","Cib2","Tph1"), 
        pt.size = 0, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2, features = c("Ccer2","Pcp4","Pou4f3","Hes6","Rprm","Psph","Atoh1","Gadd45g","Calm1","Rasd2","Tesc","Fgfr3","Hpcal1","Grp","Tmem255b","Dlk2","Cryab","Kif19a","Cib2","Tph1"),
            pt.size =  0.05, ncol = 5)


VlnPlot(E14_2, features = c("Cdkn1b","Sox2"), 
        pt.size = 0.05, ncol = 5, assay = "RNA",slot = "counts",log = T)+NoLegend()

FeaturePlot(E14_2, features = c("Cdkn1b","Sox2"),
            pt.size =  0.05, ncol = 5)

saveRDS(E14_2, file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\重复repeat\\E14 repeat\\E14_2")

setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\E14 repeat")
E14_2.named<-readRDS('E14_2.named')
new.cluster.ids_2 <- c("KO","IdC","KO","PsC","LER","OC90.Sparcl1","KO","KO","OC90/Otoa",
                       "KO","KO","HC","OC90.Sparcl1")

FeaturePlot(E14_2.named, features = c("Atoh1"),
            pt.size =  0.05, ncol = 1)


names(new.cluster.ids_2) <- levels(E14_2)
E14_2.named <- RenameIdents(E14_2, new.cluster.ids_2)
DimPlot(E14_2.named, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) + NoLegend()
saveRDS(E14_2.named, file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2\\E14 repeat\\E14_2.named")


#########select HC and SC###########
E14_select_HC_SC <- E14_2.named[,E14_2.named@meta.data$seurat_clusters %in% c(3,11)]
saveRDS(E14_select_HC_SC, file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\E14_select_HC_SC")
write.table(E14_select@assays$RNA@counts, "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2\\E14_select\\E14_select_matrix.txt",sep="\t", quote=F, col.names=NA)

#########select only SC###########
E14_select_SC <- E14_2.named[,E14_2.named@meta.data$seurat_clusters %in% c(3)]
saveRDS(E14_select_SC, file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only SC\\E14_select_SC_foldchange")
write.table(E14_select@assays$RNA@counts, "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only SC\\E14_select_only_SC_matrix.txt",sep="\t", quote=F, col.names=NA)

#########select only HC###########

E14_select_hc <- E14_2.named[,E14_2.named@meta.data$seurat_clusters %in% c(11)]

saveRDS(E14_select_hc, file = "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only HC\\E14_select_hc_foldchange")

write.table(E14_select_hc@assays$RNA@counts, "C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only HC\\E14_select_only_HC_matrix.txt",sep="\t", quote=F, col.names=NA)












#cca
E14_2_E16.cca <- SplitObject(E14_2, split.by = "time")
for (i in 1:length(E14_2_E16.cca)){
  E14_2_E16.cca[[i]]<- NormalizeData(E14_2_E16.cca[[i]],verbose = FALSE)
  E14_2_E16.cca[[i]]<- FindVariableFeatures(E14_2_E16.cca[[i]],selection.method = "vst",nfeatures = 2000, 
                                            verbose = FALSE)
}

reference.list <- E14_2_E16.cca[c("E14_2","E16_1","E16_2","E16_3")]
CCA.anchors <- FindIntegrationAnchors(object.list=reference.list,dims=1:50)

E14_2_E16.cca <-IntegrateData(anchorset=CCA.anchors,dims=1:50 )
DefaultAssay(E14_2_E16.cca)<-"integrated"
E14_2_E16.cca<-ScaleData(E14_2_E16.cca,verbose = FALSE)
E14_2_E16.cca<-RunPCA(E14_2_E16.cca,npcs = 30,verbose = FALSE)
E14_2_E16.cca<-RunUMAP(E14_2_E16.cca,reduction = "pca",dims=1:30 )
p3<-DimPlot(E14_2_E16.cca, reduction = "umap",group.by = "time")
p4<-DimPlot(E14_2_E16.cca, reduction = "umap",label = TRUE,repel = TRUE)
p3
p4
saveRDS(E14_16_integrated, file = "D:\\E16\\CCA.integrated_p60_120_E14_16")

E14_16_integrated.markers <- FindAllMarkers(E14_16_integrated, 
                                            only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
E14_16_integrated.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top30 <- E14_16_integrated.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
DoHeatmap(E14_16_integrated, features = top30$gene) + NoLegend()
cluster10.markers <- FindMarkers(E14_16_integrated, ident.1 = 10, min.pct = 0.25)
top200 <- cluster10.markers %>% top_n(n = 200, wt = avg_log2FC)
write.table(top200,"E14_16_integrated_cluster10_top200",sep = "\t")

#SC:Sox2,Cdkn1b,Prss23,Emid1 and Npy in IPCs, and Matn4 in inner phalangeal cells
VlnPlot(E14_2_E16.cca, features = c("Cdkn1b","Prss23","Sox2","Emid1","Npy","Matn4","Lfng"), 
        pt.size = 0.05, ncol = 3, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2_E16.cca, features = c("Cdkn1b","Prss23","Sox2","Emid1","Npy","Matn4","Lfng"),
            pt.size = 0.2, ncol = 3)

#HCs:Rprm, Cd164l2, Ccer2, and Gng8,Tbx2   IHC:"Pvalb","Fgf8","Otof"

VlnPlot(E14_2_E16.cca, features = c("Rprm","Cd164l2","Ccer2","Gng8","Tbx2","Pvalb","Fgf8","Otof"), 
        pt.size = 0.05, ncol = 3, assay = "RNA",slot = "counts",log = T)+NoLegend()
FeaturePlot(E14_2_E16.cca, features = c("Rprm","Cd164l2","Ccer2","Gng8","Tbx2","Pvalb","Fgf8","Otof"))


saveRDS(E14_2_E16_merge, file = "C:\\Users\\zhaoyf\\Desktop\\E14_16\\E14_2 and E16_1_2_3(seperate)\\E14_2_E16_merge")
saveRDS(E14_2, file = "C:\\Users\\zhaoyf\\Desktop\\E14_16\\E14_2 and E16_1_2_3(seperate)\\E14_2")

saveRDS(E14_2_E16.cca, file = "C:\\Users\\zhaoyf\\Desktop\\E14_16\\E14_2 and E16_1_2_3(seperate)\\E14_2_E16.cca")



