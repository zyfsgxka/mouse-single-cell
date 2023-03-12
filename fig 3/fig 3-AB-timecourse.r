#library(Mfuzz)
#library(RColorBrewer)
library(ggplot2)
library(reshape2)

########################################################################################
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\all in one")
outfolder <- ".\\2.timecourse_correlation\\"
#######################################################################################

# 1. input data #####################################################################

######## 1b. sc cell only ############
tp <- "sc_cell_P33_E14_"
## all expression matirx
dd1 <- read.csv("SC_seurat_fc_merge.txt",sep="\t")
## E14-P33 maker genes
e14.p33 <- read.csv("P33_E14_SC_integrated_fc.txt",sep="\t")
######## 1b END ######################

head(dd1)
row.names(dd1) <- dd1$gene
dd1$E14 <- 0
colnames(dd1)
d1 <- dd1[,c("E14","E16_E14","P1_E14","P7_E14","P12_E16","P33_E14")]
cname <- c("E14","E16","P1","P7","P12","P33") # group names
colnames(d1) <- cname

# 2. spearman correlation test : ####################################################
e14.p33 <- e14.p33[e14.p33$p_val_adj<0.05,]
d.cor <- d1[row.names(d1) %in% row.names(e14.p33),] # select gene subset 
d.cor <- d1[abs(d1$P33)>1,] # select gene subset 

time <- c(0.14,0.16,1,7,12,33)  # dependent factor for correlation (development stage)
#time <- c(1,2,3,4,5,6)

## calculate spearman correlation coefficient ##
  p.sp <- c() 
  r.sp <- c()  
  for(i in 1:nrow(d.cor)){
    sp <- cor.test(x=time,y=as.numeric(d.cor[i,]),method = "spearman")
    p.sp <- c(p.sp,sp$p.value)
    r.sp <- c(r.sp,as.numeric(sp$estimate))
  }
  d.cor$r.sp <- r.sp   # correlation coefficient
  d.cor$p.sp <- p.sp   # p value

## candiate gene selection ##
dd <- d.cor[abs(d.cor$r.sp)>0.9 & d.cor$p.sp<0.01 ,]  # candidate genes
heatmap(as.matrix(dd[,1:6]),Colv = NA,scale = "none",labCol = cname)

## expression data normalization: ##
head(dd)
ddm <- as.matrix(dd[,1:6])
mtx <- t(scale(t(ddm)))  # center by row
heatmap(mtx,Colv = NA,scale = "none",labCol = cname)

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)


Heatmap(mtx,col=colorRampPalette(rev(c("#Ff8c00","#b0e0e6")))(102),show_row_names = FALSE,column_order = NULL)


#library(pheatmap)
#pheatmap(mtx, cluster_cols  = FALSE, cluster_rows = TRUE,Colv = NA,scale = "none",
         #color=colorRampPalette(rev(c("red","yellow")))(256),border_color = F,show_row_names = FALSE)


#"#ffd700","#f4a460","#b0e0e6"
                                                                                                                                                                                  
k.num <- 2  ## set number of  clusters 

k.data <- mtx 
kmeans_clustering <- kmeans(k.data, centers=k.num, iter.max=100, nstart=5)
clust <- as.data.frame(kmeans_clustering$cluster)
kmeans.clusterd <- cbind(as.data.frame(k.data),clust)
d.cluster <- cbind(dd,clust)
dd$cluster <- clust$`kmeans_clustering$cluster`

## plot clusters ##
#### a. line plot ###
plots_per_row <- 2
plots_per_col <- round(k.num/plots_per_row)+1

#pdf(file = paste(outfolder,tp,nrow(dd),"kmeans_",k.num,"cluster_line_plot.pdf",sep=""),height = round(k.num/3+0.4)*3)
#dev.off()
par(mfrow=c(plots_per_col, plots_per_row))
par(cex=0.6)
par(mar=c(2,4,4,2))
line.col <- rgb(221,221,221,alpha = 150,maxColorValue = 255) # backgournd line color
mean.col <- "coral"  # mean line color
line.col <- "#dcdcdc"
mean.col <- "#Ff8c00"  # mean line color


ylab <- "RNA level"

for (i in 1:k.num) {
  data.ksub <-  kmeans.clusterd[kmeans.clusterd[,length(colnames(kmeans.clusterd))]==i,]
  data <- data.ksub[,1:length(colnames(data.ksub))-1]
  ymin <-  min(data); ymax <-  max(data)
  plot_label <-  paste("cluster_", i,"_", length(data[,1]), "gene", sep='')
  
  # basic plot frame
  plot(as.numeric(data[1,]), type='l', ylim=c(ymin,ymax), main=plot_label, 
       col=line.col, xaxt='n', xlab='',ylab=ylab,lwd=0.5,)
  axis(side=1, at=1:length(data[1,]), labels=cname, las=2)
  # add more lines
  for(r in 2:length(data[,1])) {
    points(as.numeric(data[r,]), type='l', col=line.col )
  }
  # add means
  points(as.numeric(colMeans(data)), type='o', col=mean.col,lwd=3) 
}


write.table(dd,file = paste(outfolder,tp,nrow(dd),"kmeans_",k.num,"cluster.txt",sep=""),quote = F,sep = "\t")
for(i in 1:k.num){
  write.table(dd[dd$cluster==i,],file = paste(outfolder,tp,nrow(dd),"kmeans_",k.num,"cluster",i,".txt",sep=""),quote = F,sep = "\t")
}



