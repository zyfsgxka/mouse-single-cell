#########################################################################################
###############################Umap and featureplot######################################
##########################################################################################

library(Seurat)
#######################E14 repeat#########################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\all in one")
E14_2.named<-readRDS('E14_2.named')
new.cluster.ids_2 <- c("KO","IdC","KO","PsC","LER","OC90.Sparcl1","KO","KO","OC90/Otoa",
                       "KO","KO","HC","OC90.Sparcl1")
#FeaturePlot(E14_2.named, features = c("Atoh1"),pt.size =  0.05, ncol = 1)
#names(new.cluster.ids_2) <- levels(E14_2)
E14_2.named <- RenameIdents(E14_2.named, new.cluster.ids_2)
DimPlot(E14_2.named, reduction = "umap", label = TRUE, pt.size = 0.5,label.size =7) +ggtitle("E14")+ NoLegend()
#####umap Plot export size:8x6
FeaturePlot(E14_2.named, features = c("Atoh1"),
            pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(E14_2.named, features = c("Sox2"),
            pt.size = 1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(E14_2.named, features = c("Pcdh15"),
            pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(E14_2.named, features = c("Ezh2"),
            pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))

FeaturePlot(E14_2.named, features = c("Klf15"),
            pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
#####FeaturePlot export size:4x5



########################E16 repeat######################

setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\all in one\\umap and feature")
E16_cca<-readRDS("E16_cca")
#E16_cca.markers <- FindAllMarkers(E16_cca, min.pct = 0.25, logfc.threshold = 0.25)
#write.table(E16_cca.markers,"E16_cca.all_markers.xls",sep = "\t",quote = F,row.names = F)
new.cluster.ids_3 <- c("KO","ISC","L.KO","PsC","Igfbp2+","IdC","LER",
                       "HC","OC90","Ube2c+","Malat1+","Ube2c+","Aspa+","IPhC",
                       "Gm43376+","Hensen","IPC","Trpm1+",
                       "CC/OSC","CC/OSC")
names(new.cluster.ids_3) <- levels(E16_cca)
cca.renamed <- RenameIdents(E16_cca, new.cluster.ids_3)
p5<-DimPlot(cca.renamed, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 7) +ggtitle("E16")+NoLegend()
p5
FeaturePlot(cca.renamed, features = c("Atoh1"),pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(cca.renamed, features = c("Sox2"),
            pt.size = 1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(cca.renamed, features = c("Pcdh15"),
            pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))




###############P1 repeat################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\all in one\\umap and feature")
P1_cca<-readRDS("P1.cca")
DimPlot(P1_cca, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 7) +ggtitle("P1_named")+NoLegend()

#P1_cca.markers <- FindAllMarkers(P1_cca, min.pct = 0.25, logfc.threshold = 0.25)

#write.table(P1_cca.markers,"P1_cca.all_markers.xls",sep = "\t",quote = F,row.names = F)

new.cluster.ids_3 <- c("M.KO","L.KO","L.KO","Sptssa+","Mgp+","ISC/IdC","Otor+",
                       "Otor+","Ccn3+","HC","DC/OPC","Arpc1b+","CC/OSC/HeC",
                       "IPhc","Adgre1+","Hmgb2+",
                       "OC90","S100a1+","Gm14635+","HC","IPC","Trpm1+","Cpne4+",
                       "Esrrb+")
names(new.cluster.ids_3) <- levels(P1_cca)
cca.renamed <- RenameIdents(P1_cca, new.cluster.ids_3)
p5<-DimPlot(cca.renamed, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 7) +ggtitle("P1")+NoLegend()
p5
FeaturePlot(cca.renamed, features = c("Atoh1"),pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(cca.renamed, features = c("Sox2"),
            pt.size = 1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(cca.renamed, features = c("Pcdh15"),
            pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))



###########################P7 repeat##################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\all in one\\umap and feature")
P7_cca<-readRDS("P7_cca")
#P7_cca.markers <- FindAllMarkers(P7_cca, min.pct = 0.25, logfc.threshold = 0.25)
#write.table(P7_cca.markers,"P7_cca.all_markers.xls",sep = "\t",quote = F,row.names = F)
new.cluster.ids_3 <- c("L.KO","Fbxo2+","Apod+","M.KO","Ccn3+","Rarres1+","IPhC","CC/OS","HeC/Glia",
                       "OPC/IPC/DC","HC","C1qa+","L.KO","HC","Aspn+")
names(new.cluster.ids_3) <- levels(P7_cca)
cca.renamed <- RenameIdents(P7_cca, new.cluster.ids_3)
p5<-DimPlot(cca.renamed, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 7) +ggtitle("P7")+NoLegend()
p5
FeaturePlot(cca.renamed, features = c("Atoh1"),pt.size =  2, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(cca.renamed, features = c("Sox2"),
            pt.size = 1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(cca.renamed, features = c("Pcdh15"),
            pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))




##############P12 repeat###################
cca.object<-readRDS("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\all in one\\umap and feature\\P12_cca")
new.cluster.ids_2 <- c("Notum+","SC","Slfn14+","Sdpr+","Itgam+","Opcml+","HC",
                       "HC","HC","SC","Bmp8a+","Mxd3+")
names(new.cluster.ids_2) <- levels(cca.object)
cca.named <- RenameIdents(cca.object, new.cluster.ids_2)
DimPlot(cca.named, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 7) + ggtitle("P12")+NoLegend()
FeaturePlot(cca.named, features = c("Atoh1"),pt.size =  2, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(cca.named, features = c("Sox2"),
            pt.size = 1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(cca.named, features = c("Pcdh15"),
            pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))



#################P33 repeat##################3
cca.named<-readRDS("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\all in one\\P33_cca.renamed")
s.object<-cca.named
new.cluster.ids_2 <- c("Alas2+","HC","SC","SC","Cldn3+","SC","HC")

names(new.cluster.ids_2) <- levels(s.object)
cca.named <- RenameIdents(s.object, new.cluster.ids_2)
DimPlot(cca.named, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 7) +ggtitle("P33")+ NoLegend()
FeaturePlot(cca.named, features = c("Ezh2"),pt.size =  2, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(cca.named, features = c("Sox2"),
            pt.size = 1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))
FeaturePlot(cca.named, features = c("Pcdh15"),
            pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))


FeaturePlot(cca.named, features = c("Ezh2"),
           pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))


FeaturePlot(cca.named, features = c("Klf15"),
            pt.size =  1, ncol = 1,cols=c("#dcdcdc","#Ff8c00","#ff4500"))


#######################################################################
#############Atoh1,sox2,pcdh15 changing curve#########################
########################################################################
library(reshape2)     
library(ggplot2)
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\all in one\\SC foldchange noslots")
all_FC_merge <-  read.csv("SC_seurat_fc_merge.txt", header=T, sep = "\t")
all_FC_merge$E14 <- 0

#all_FC_merge$gene <-rownames(all_FC_merge)
colnames(all_FC_merge)<-c("gene","E16","P1","P7","P12","P33","E14")
all_FC_merge <- all_FC_merge[,c(1,7,2,3,4,5,6)]
####################################
target_gene<-c("Apoe")
all_FC_merge<- all_FC_merge[all_FC_merge$gene%in%target_gene,c(1,2,3,4,5,6,7)] ##get names of test pairs   
aql <- melt(all_FC_merge)

name<-paste(aql$gene,sep = "")

ggplot(aql,aes(variable,value,group =gene,color=gene)) +geom_line(size=1.5)+ylab("logFC")+theme_bw() +facet_wrap(~gene)+
  theme_classic()+theme(legend.position = 'none')+
  geom_point(size=5, shape=20)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) +ylim(-5,5)+
  scale_colour_manual(values=c("#Ff8c00"))+  
  theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),
  axis.title.x= element_text(size=15),axis.title.y= element_text(size=15),
  text = element_text(size=30) )

#######################################
target_gene<-c("Ezh2","Klf15")
all_FC_merge<- all_FC_merge[all_FC_merge$gene%in%target_gene,c(1,2,3,4,5,6,7)] ##get names of test pairs   
aql <- melt(all_FC_merge)

name<-paste(aql$gene,sep = "")

ggplot(aql,aes(variable,value,group=gene,color=gene)) +geom_line(size=1.5)+ylab("logFC") +
  theme_classic()+theme(legend.position = 'right')+
  geom_point(size=5, shape=20)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) +ylim(-0.7,0.7)+
  scale_colour_manual(values=c("#b0e0e6","#Ff8c00"))+  
  theme(axis.text.x = element_text(size=15),axis.text.y=element_text(size=15),
        axis.title.x= element_text(size=15),axis.title.y= element_text(size=15),
        text = element_text(size=15) )
############export size:6x5













##################################################################################
#################################cell  propotion ################################
#################################################################################

library(Seurat)
library(ggplot2)
library(ggsci)

##################E14####################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\E14 repeat")
E14_2.named<-readRDS('E14_2')
new.cluster.ids_3 <- c("other cells","other cells","other cells","PsC","other cells","other cells","other cells","other cells","other cells",
                       "other cells","other cells","HC","other cells")
names(new.cluster.ids_3) <- levels(E14_2.named)
E14_2.named <- RenameIdents(E14_2.named, new.cluster.ids_3)
E14_select <- E14_2.named[,E14_2.named@meta.data$seurat_clusters %in% c(3,11)]

table(Idents(E14_select))###Count the number of cells in each cluster1
prop.table(table(Idents(E14_select)))###Calculating Cell Proportion
cell.prop_E14<-as.data.frame(prop.table(table(Idents(E14_select))))
#cell.prop_E14<-as.data.frame(table(Idents(E14)))
colnames(cell.prop_E14)<-c("cell_type","proportion")
cell.prop_E14$time<-"E14"
signif(cell.prop_E14$proportion, 4)######Keep four decimals
ggplot(cell.prop_E14,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_igv()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.02), position=position_dodge(1), vjust=0)+xlab("E14")


##################E16####################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\E16 repeat")

E16_cca<-readRDS("E16_cca")


new.cluster.ids_3 <- c("other cells","other cells","other cells","PsC","other cells","other cells","other cells",
                       "HC","other cells","other cells","other cells","other cells","other cells","SC",
                       "other cells","SC","SC","other cells",
                       "SC","SC")
names(new.cluster.ids_3) <- levels(E16_cca)
E16 <- RenameIdents(E16_cca, new.cluster.ids_3)
E16_select <- E16[,E16@meta.data$seurat_clusters %in% c(3,7,13,15,16,18,19)]
#
table(Idents(E16_select))
prop.table(table(Idents(E16_select)))
cell.prop_E16<-as.data.frame(prop.table(table(Idents(E16_select))))
#cell.prop_E16<-as.data.frame(table(Idents(E16)))
colnames(cell.prop_E16)<-c("cell_type","proportion")
cell.prop_E16$time<-"E16"
signif(cell.prop_E16$proportion, 4)

#
ggplot(cell.prop_E16,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_igv()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.005), position=position_dodge(1), vjust=0)

############### P1  #########


setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P1_repeat")
P1<-readRDS('P1.cca')

new.cluster.ids_3 <- c("other cells","other cells","other cells","other cells","other cells","other cells","other cells",
                       "other cells","other cells","HC","SC","other cells","SC",
                       "SC","other cells","other cells",
                       "other cells","other cells","other cells","HC","SC","other cells","other cells",
                       "other cells")
names(new.cluster.ids_3) <- levels(P1)
P1 <- RenameIdents(P1, new.cluster.ids_3)
P1_select <- P1[,P1@meta.data$seurat_clusters %in% c(9,10,12,13,19,20)]
#
table(Idents(P1_select))
prop.table(table(Idents(P1_select)))
cell.prop_P1<-as.data.frame(prop.table(table(Idents(P1_select))))
#cell.prop_P1<-as.data.frame(table(Idents(P1)))
colnames(cell.prop_P1)<-c("cell_type","proportion")
cell.prop_P1$time<-"P1"
signif(cell.prop_P1$proportion, 4)

ggplot(cell.prop_P1,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_igv()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.005), position=position_dodge(1), vjust=0)

############ P7
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P7_repeat")
P7_cca<-readRDS("P7_cca")
new.cluster.ids_3 <- c("other cells","other cells","other cells","other cells","other cells","other cells","SC","SC","SC",
                       "SC","HC","other cells","other cells","HC","other cells")
names(new.cluster.ids_3) <- levels(P7_cca)
P7 <- RenameIdents(P7_cca, new.cluster.ids_3)
P7_select <- P7[,P7@meta.data$seurat_clusters %in% c(6,7,8,9,10,13)]
#
table(Idents(P7_select))
prop.table(table(Idents(P7_select)))
cell.prop_P7<-as.data.frame(prop.table(table(Idents(P7_select))))
#cell.prop_P7<-as.data.frame(table(Idents(P7)))
colnames(cell.prop_P7)<-c("cell_type","proportion")
cell.prop_P7$time<-"P7"
signif(cell.prop_P7$proportion, 4)

ggplot(cell.prop_P7,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_igv()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.005), position=position_dodge(1), vjust=0)


############ P12
P12<-readRDS("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P12\\P12.cca")
new.cluster.ids_2 <- c("other cells","SC","other cells","other cells","other cells","other cells","HC",
                       "HC","HC","SC","other cells","other cells")
names(new.cluster.ids_2) <- levels(P12)
P12 <- RenameIdents(P12, new.cluster.ids_2)
P12_select <- P12[,P12@meta.data$seurat_clusters %in% c(1,6,7,8,9)]
#
table(Idents(P12_select))
#prop.table(table(Idents(P12)))
cell.prop_P12<-as.data.frame(prop.table(table(Idents(P12_select))))
#cell.prop_P12<-as.data.frame(table(Idents(P12)))
colnames(cell.prop_P12)<-c("cell_type","proportion")
cell.prop_P12$time<-"P12"
signif(cell.prop_P12$proportion, 4)

ggplot(cell.prop_P12,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_igv()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.005), position=position_dodge(1), vjust=0)


############ P33
P33<-readRDS("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P33\\P33.object")

new.cluster.ids_2 <- c("other cells","HC","SC","SC","Cldn3+","SC","HC")
names(new.cluster.ids_2) <- levels(P33)
P33 <- RenameIdents(P33, new.cluster.ids_2)
P33_select <- P33[,P33@meta.data$seurat_clusters %in% c(1,2,3,5,6)]
#
table(Idents(P33_select))
#prop.table(table(Idents(P33)))
cell.prop_P33<-as.data.frame(prop.table(table(Idents(P33_select))))
#cell.prop_P33<-as.data.frame(table(Idents(P33)))
colnames(cell.prop_P33)<-c("cell_type","proportion")
cell.prop_P33$time<-"P33"
signif(cell.prop_P33$proportion, 4)

ggplot(cell.prop_P33,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_igv()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.005), position=position_dodge(1), vjust=0)


######### merge#############
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\all in one\\cell proportion")
merge<- rbind(cell.prop_E14,cell.prop_E16,cell.prop_P1,cell.prop_P7,cell.prop_P12,cell.prop_P33)
merge$time<-factor(merge$time,levels=c("E14","E16","P1","P7","P12","P33"),ordered=TRUE)
merge$cell_type<-factor(merge$cell_type,levels=c("PsC","SC","HC"),ordered=TRUE)

ggplot(merge,aes(time,proportion,fill=cell_type))+geom_bar(stat='identity',position='stack',width = 0.9)+theme_classic()+
  ggtitle("")+theme_bw()+#theme(axis.ticks.length=unit(0.5,'cm'))
  scale_fill_manual(values=c("#ffd700","#f4a460","#b0e0e6"))+
  theme(panel.grid =element_blank())
  #theme(axis.text.x = element_text(size=18),axis.text.y=element_text(size=18),
     # axis.title.x= element_text(size=18),axis.title.y= element_text(size=18),
      #text = element_text(size=18) )











#############################################################################################
################################### GSEA######################################
##############################################################################################
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(KEGG.db)
library(enrichplot)
library(dplyr)
##########    set input file and threshold accordingly:   #############################################


time <- "cochlearP7-E14_all_exact_"

#####################################################################################################



setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\all in one\\umap and feature\\P33-GSEA")
#all_id<-read.table("all_id.txt",sep = "\t",header = T)
#all_id <- na.omit(all_id)
dd <- read.table("P33_E14_SC_integrated_fc.txt",header=T,sep = "\t")
dd$gene<-rownames(dd)
colnames(dd)
#colnames(dd)<-c("logFC","logCPM","PValue","FDR","gene")


keytypes(org.Mm.eg.db)
ENTREZID <- bitr(dd$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
ENTREZID$ENTREZID
colnames(ENTREZID)<-c("gene","ENTREZID")


d2 <- merge(dd,ENTREZID,by.x="gene")
d2<- arrange(d2,desc(avg_log2FC))

## 1. GSEA: KEGG  #######################################################
# rank genelist
d2$fcsign <- sign(d2$avg_log2FC)
d2$logP <- -log10(d2$p_val_adj)
d2$metric <-  d2$logP/d2$fcsign
GSEA_input<-d2$avg_log2FC
# GSEA_input<-d2$logFC
names(GSEA_input) = as.character(d2$ENTREZID)
# GSEA_input<-d2[!is.na(d2$ENTREZID),]$logFC
# names(GSEA_input) = as.character(d2[!is.na(d2$ENTREZID),]$ENTREZID)
GSEA_input = sort(GSEA_input, decreasing = T)
GSEA_input
gseKEGG.res <- gseKEGG(GSEA_input, organism = "mmu",keyType = "ncbi-geneid",nPerm = 1000, minGSSize = 5, maxGSSize = 1000, pvalueCutoff=1)
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only SC\\foldchange\\P7-E14 gsea")
write.table(as.data.frame(gseKEGG.res),file=paste(time,"GSEA.kegg.enrich.ncbi.txt",sep = "_"),quote = F,sep="\t",row.names = F)

inpath<-c("mmu04310","mmu04064","mmu04115","mmu04010","mmu04330","mmu04151","mmu04668",
          "mmu00190","mmu04210","mmu00010","mmu04620","mmu04621","mmu04514")

inpath1 <-"mmu04310"	#Wnt signaling pathway
inpath2<-"mmu04064" #	NF-kappa B signaling pathway   
inpath3<-"mmu04115" # 	p53 signaling pathway
inpath4<-"mmu04010" #	MAPK signaling pathway
inpath5<-"mmu04330" #	Notch signaling pathway
inpath6<-"mmu04151" #	PI3K-Akt signaling pathway
inpath7<-"mmu04668" #	TNF signaling pathway
mmu00190
inpath8<-"mmu00190" #	Oxidative phosphorylation
inpath9<-"mmu04210" #	Apoptosis
inpath10<-"mmu00010" #	Glycolysis / Gluconeogenesis
inpath11<-"mmu04620"#Toll-like receptor signaling pathway
inpath12<-"mmu04621"#NOD-like receptor signaling pathway
inpath13<-"mmu04514"#cell adhesion

ridgeplot(inpath5)

for (i in 1:length(inpath)){
  inpath_w<-inpath[i]
  p1 <- gseaplot2(gseKEGG.res,geneSetID = inpath5,pvalue_table = F,color = "#f4a460",title = "Notch signaling pathway",
                  base_size = 15,subplots=c(1,2),rel_heights = c(2,1)) 
  p1    
  ggsave(paste(time,inpath_w,"GSEA_kegg.pdf",sep="_"),p1,width = 10,height = 5) 
}

#core <- strsplit(gseKEGG.res@result[gseKEGG.res@result$ID==inpath2,]$core_enrichment,split = "/")[[1]]

#bitr(geneID = core,fromType = "ENTREZID",toType ="SYMBOL" ,OrgDb = org.Mm.eg.db)


##### GSEA: GO

gseGO.res <- gseGO(GSEA_input,ont = "BP",OrgDb = org.Mm.eg.db,keyType = "ENTREZID",pvalueCutoff = 1 )
write.table(as.data.frame(gseGO.res),file=paste(time,"GSEA.GOBP.enrich.ncbi.txt",sep = "_"),quote = F,sep="\t",row.names = F)

inpath<-c("GO:1901796","GO:0030177","GO:0043122","GO:0043410","GO:0007219")
inpath1 <-"GO:1901796" #	regulation of signal transduction by p53 class mediator
inpath2 <-"GO:0030177" #		positive regulation of Wnt signaling pathway
inpath3 <-"GO:0043122"	#regulation of I-kappaB kinase/NF-kappaB signaling
inpath4 <-"GO:0043410"	# 	positive regulation of MAPK cascade
inpath5 <-"GO:0007219" #		Notch signaling pathway
for (i in 1:length(inpath)){
  inpath_w2<-inpath[i]
  p2 <- gseaplot2(gseGO.res,geneSetID = inpath_w2,pvalue_table = T)
  p2
  ggsave(paste(time,strsplit(inpath_w2,":")[[1]][2],"GSEA.GOBP.pdf",sep="_"),p2,width =10,height = 5)
}












#################################################################################
#################################cell proportion auditory cell #####################################
##################################################################################

library(Seurat)
library(ggplot2)
library(ggsci)

##################E14####################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\E14 repeat")
E14_2.named<-readRDS('E14_2')
#new.cluster.ids_3 <- c("other cells","other cells","other cells","PsC","other cells","other cells","other cells","other cells","other cells",
"other cells","other cells","HC","other cells")
new.cluster.ids_3 <- c("other cells","other cells","other cells","sensory cells",
                       "other cells","other cells","other cells","other cells","other cells",
                       "other cells","other cells","sensory cells","other cells")
names(new.cluster.ids_3) <- levels(E14_2.named)
E14_2.named <- RenameIdents(E14_2.named, new.cluster.ids_3)
E14_select <- E14_2.named#[,E14_2.named@meta.data$seurat_clusters %in% c(3,11)]

table(Idents(E14_select))###Count the number of cells in each cluster1
prop.table(table(Idents(E14_select)))###Calculating Cell Proportion
cell.prop_E14<-as.data.frame(prop.table(table(Idents(E14_select))))
#cell.prop_E14<-as.data.frame(table(Idents(E14)))
colnames(cell.prop_E14)<-c("cell_type","proportion")
cell.prop_E14$time<-"E14"
signif(cell.prop_E14$proportion, 4)######Keep four decimals
ggplot(cell.prop_E14,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_d3()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.02), position=position_dodge(1), vjust=0)+xlab("E14")


##################E16####################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\E16 repeat")

E16_cca<-readRDS("E16_cca")


#new.cluster.ids_3 <- c("other cells","other cells","other cells","PsC","other cells","other cells","other cells",
"HC","other cells","other cells","other cells","other cells","other cells","SC",
"other cells","SC","SC","other cells",
"SC","SC")
new.cluster.ids_3 <- c("other cells","other cells","other cells","sensory cells",
                       "other cells","other cells","other cells",
                       "sensory cells","other cells","other cells","other cells",
                       "other cells","other cells","sensory cells",
                       "other cells","sensory cells","sensory cells","other cells",
                       "sensory cells","sensory cells")
names(new.cluster.ids_3) <- levels(E16_cca)
E16 <- RenameIdents(E16_cca, new.cluster.ids_3)
E16_select <- E16#[,E16@meta.data$seurat_clusters %in% c(3,7,13,15,16,18,19)]
#
table(Idents(E16_select))
prop.table(table(Idents(E16_select)))
cell.prop_E16<-as.data.frame(prop.table(table(Idents(E16_select))))
#cell.prop_E16<-as.data.frame(table(Idents(E16)))
colnames(cell.prop_E16)<-c("cell_type","proportion")
cell.prop_E16$time<-"E16"
signif(cell.prop_E16$proportion, 4)

#
ggplot(cell.prop_E16,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_igv()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.005), position=position_dodge(1), vjust=0)

############### P1  #########
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P1_repeat")
P1<-readRDS('P1.cca')

#new.cluster.ids_3 <- c("other cells","other cells","other cells","other cells","other cells","other cells","other cells",
"other cells","other cells","HC","SC","other cells","SC",
"SC","other cells","other cells",
"other cells","other cells","other cells","HC","SC","other cells","other cells",
"other cells")
new.cluster.ids_3 <- c("other cells","other cells","other cells",
                       "other cells","other cells","other cells","other cells",
                       "other cells","other cells","sensory cells","sensory cells","other cells","sensory cells",
                       "sensory cells","other cells","other cells",
                       "other cells","other cells","other cells","sensory cells","sensory cells","other cells","other cells",
                       "other cells")
names(new.cluster.ids_3) <- levels(P1)
P1 <- RenameIdents(P1, new.cluster.ids_3)
P1_select <- P1#[,P1@meta.data$seurat_clusters %in% c(9,10,12,13,19,20)]
#
table(Idents(P1_select))
prop.table(table(Idents(P1_select)))
cell.prop_P1<-as.data.frame(prop.table(table(Idents(P1_select))))
#cell.prop_P1<-as.data.frame(table(Idents(P1)))
colnames(cell.prop_P1)<-c("cell_type","proportion")
cell.prop_P1$time<-"P1"
signif(cell.prop_P1$proportion, 4)

ggplot(cell.prop_P1,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_igv()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.005), position=position_dodge(1), vjust=0)

############ P7
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P7_repeat")
P7_cca<-readRDS("P7_cca")
#new.cluster.ids_3 <- c("other cells","other cells","other cells","other cells","other cells","other cells","SC","SC","SC",
"SC","HC","other cells","other cells","HC","other cells")
new.cluster.ids_3 <- c("other cells","other cells","other cells","other cells",
                       "other cells","other cells","sensory cells","sensory cells","sensory cells",
                       "sensory cells","sensory cells","other cells","other cells","sensory cells","other cells")
names(new.cluster.ids_3) <- levels(P7_cca)
P7 <- RenameIdents(P7_cca, new.cluster.ids_3)
P7_select <- P7#[,P7@meta.data$seurat_clusters %in% c(6,7,8,9,10,13)]
#
table(Idents(P7_select))
prop.table(table(Idents(P7_select)))
cell.prop_P7<-as.data.frame(prop.table(table(Idents(P7_select))))
#cell.prop_P7<-as.data.frame(table(Idents(P7)))
colnames(cell.prop_P7)<-c("cell_type","proportion")
cell.prop_P7$time<-"P7"
signif(cell.prop_P7$proportion, 4)

ggplot(cell.prop_P7,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_igv()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.005), position=position_dodge(1), vjust=0)


############ P12
P12<-readRDS("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P12\\P12.cca")
#new.cluster.ids_2 <- c("other cells","SC","other cells","other cells","other cells","other cells","HC",
"HC","HC","SC","other cells","other cells")
new.cluster.ids_2 <- c("other cells","sensory cells","other cells","other cells",
                       "other cells","other cells","sensory cells",
                       "sensory cells","sensory cells","sensory cells","other cells","other cells")
names(new.cluster.ids_2) <- levels(P12)
P12 <- RenameIdents(P12, new.cluster.ids_2)
P12_select <- P12#[,P12@meta.data$seurat_clusters %in% c(1,6,7,8,9)]
#
table(Idents(P12_select))
#prop.table(table(Idents(P12)))
cell.prop_P12<-as.data.frame(prop.table(table(Idents(P12_select))))
#cell.prop_P12<-as.data.frame(table(Idents(P12)))
colnames(cell.prop_P12)<-c("cell_type","proportion")
cell.prop_P12$time<-"P12"
signif(cell.prop_P12$proportion, 4)

ggplot(cell.prop_P12,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_igv()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.005), position=position_dodge(1), vjust=0)


############ P33
P33<-readRDS("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\P33\\P33.object")

#new.cluster.ids_2 <- c("other cells","HC","SC","SC","other cells","SC","HC")
new.cluster.ids_2 <- c("other cells","sensory cells","sensory cells","sensory cells",
                       "other cells","sensory cells","sensory cells")
names(new.cluster.ids_2) <- levels(P33)

P33 <- RenameIdents(P33, new.cluster.ids_2)
P33_select <- P33#[,P33@meta.data$seurat_clusters %in% c(1,2,3,5,6)]
#
table(Idents(P33_select))
#prop.table(table(Idents(P33)))
cell.prop_P33<-as.data.frame(prop.table(table(Idents(P33_select))))
#cell.prop_P33<-as.data.frame(table(Idents(P33)))
colnames(cell.prop_P33)<-c("cell_type","proportion")
cell.prop_P33$time<-"P33"
signif(cell.prop_P33$proportion, 4)

ggplot(cell.prop_P33,aes(cell_type,proportion,fill=cell_type))+geom_bar(stat="identity",position = position_dodge(1))+
  ggtitle("")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+scale_fill_igv()+
  geom_text(aes(label=signif(proportion, 4), y=proportion+0.005), position=position_dodge(1), vjust=0)


######### merge
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\all in one\\cell proportion")
merge<- rbind(cell.prop_E14,cell.prop_E16,cell.prop_P1,cell.prop_P7,cell.prop_P12,cell.prop_P33)
merge$time<-factor(merge$time,levels=c("E14","E16","P1","P7","P12","P33"),ordered=TRUE)
#merge$cell_type<-factor(merge$cell_type,levels=c("other cells","PsC","SC","HC"),ordered=TRUE)
merge$cell_type<-factor(merge$cell_type,levels=c("other cells","sensory cells"),ordered=TRUE)
ggplot(merge,aes(time,proportion,fill=cell_type))+geom_bar(stat='identity',position='stack',width = 0.9)+theme_classic()+
  ggtitle("")+theme_bw()+#+theme(axis.ticks.length=unit(0.5,'cm')
  scale_fill_manual(values=c("#dcdcdc","#Ff8c00"))+
  theme(panel.grid =element_blank())













