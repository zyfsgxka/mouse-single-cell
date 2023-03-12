
#############Notch target gene heatmap#############

setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\repeat2 mouse\\select only SC\\foldchange\\SC foldchange noslots")
SC_seurat_fc_merge <- read.table("SC_seurat_fc_merge.txt",header=T,sep = "\t")

target <- c("Casp3","Ccnd1","Ccne1","Cdkn1a","Cdkn1a","Ctnnb1","Dll4","Dtx1","Efnb2","Gata2","Hes5","Hey1","Hey2","Nfkb1","Nrarp"
            ,"Pax7","Rbpj","Runx2","Scgb3a2","Smarcc1","Sox9","Sp7","Spi1","Tcf3","Gzmb","Hes1","Notch3","Foxp3","Ptcra","Dlk1","Hes1")

x <- SC_seurat_fc_merge[SC_seurat_fc_merge$gene %in% target,]

head(x)
row.names(x) <- x$gene
x$E14 <- 0
colnames(x)
d1 <- x[,c("E14","E16_E14","P1_E14","P7_E14","P12_E16","P33_E14")]
cname <- c("E14","E16","P1","P7","P12","P33") # group names
colnames(d1) <- cname
mtx <- t(scale(t(d1)))
library(ComplexHeatmap)
Heatmap(mtx,col=colorRampPalette(rev(c("#6495ED","#E0FFFF")))(102),show_row_names = T,
        column_order = NULL, cluster_rows = T,name = "log2FC")

