
#############Notch GSEA core gene heatmap#############

setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only SC\\foldchange\\SC foldchange noslots")
SC_seurat_fc_merge <- read.table("SC_seurat_fc_merge.txt",header=T,sep = "\t")

############intersection############
target1 <- c("Ctbp1","Ctbp2","Kat2a","Hdac2","Hes1","Hey1","Hey2","Jag1","Notch1","Tle1","Aph1a","Hdac1","Snw1")

#############P33-P7 specific########
target2 <- c("Adam17","Numbl","Psen1","Tle3","Heyl")


#############P7-E14 specific########
target3 <- c("Crebbp","Dvl3","Lfng","Mfng","Notch3","Rbpj","Atxn1","Tle4","Maml2","Maml3","Psenen","Cir1","Dtx2")


###choose target1,2,3 as object######
x <- SC_seurat_fc_merge[SC_seurat_fc_merge$gene %in% target1,]  

head(x)
row.names(x) <- x$gene
x$E14 <- 0
colnames(x)
d1 <- x[,c("E14","E16_E14","P1_E14","P7_E14","P12_E16","P33_E14")]
cname <- c("E14","E16","P1","P7","P12","P33") # group names
colnames(d1) <- cname
mtx <- t(scale(t(d1)))
library(ComplexHeatmap)
Heatmap(mtx,col=colorRampPalette(rev(c("#FFA500","#FAF0E6")))(102),show_row_names = T,
        column_order = NULL, cluster_rows = T,name = "log2FC")




