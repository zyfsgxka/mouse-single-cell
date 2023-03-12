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
