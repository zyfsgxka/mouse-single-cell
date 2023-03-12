library(ggplot2)

setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面")
kegg <-read.table("GO target2.6.txt",sep = "\t",header = T)

colnames(kegg)

#kegg$GO_biological_process <- sapply(strsplit(as.character(kegg$GO.biological.process.complete),"(",fixed= T), "[", 1)

kegg$GO_biological_process <- sapply(strsplit(as.character(kegg$GO),"(",fixed= T), "[", 1)

colnames(kegg) <- c("GO.biological.process.complete","REFLIST","Count",
                    "expected","under","Fold.Enrichment","P-value","FDR","GO_biological_process")

kegg$GO_biological_process<-factor(kegg$GO_biological_process,levels=c("positive regulation of epithelial cell migration",
                                                                       "regulation of epithelial cell proliferation", 
                                                                       "positive regulation of cell population proliferation", "sensory organ development",
                                                                       "regulation of cell differentiation",  
                                                                       "regulation of developmental process", 
                                                                       "regulation of immune system process",   
                                                                       "positive regulation of cellular component organization",
                                                                       "nervous system development","animal organ development"),ordered=TRUE)

ggplot(kegg,aes(Fold.Enrichment,GO_biological_process))+ 
  geom_point(aes(size=Count))+
  geom_point(aes(size=Count,color=-1*log10(FDR)))+
  scale_color_gradient(low="#b0e0e6",high = "#Ff8c00")+
  labs(x="fold Enrichment",y="GO biological process",
       size="Gene number",color="-log10(FDR)")+theme_classic()  +
  theme(axis.text.x = element_text(size=13),axis.text.y=element_text(size=13),
        axis.title.x= element_text(size=13),axis.title.y= element_text(size=13))#scale_fill_manual(values=c("#478bA2","#DDC1C6","#F2A490","#E9765B"))

ggplot(kegg,aes(-1*log10(FDR),GO_biological_process))+
  geom_bar(stat="identity",fill = "#CC79A7")+ 
  labs(x="-log10(FDR)",y="GO biological process")+
  theme_classic()+geom_text(aes(label=Count),hjust = -1)
