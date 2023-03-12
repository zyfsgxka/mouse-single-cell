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