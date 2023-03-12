library(RcisTarget)


########## input data #######################################################################################
# gene list file:
setwd("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\Rcistarget_database")
d <- read.csv("C:\\Users\\zhaoyf\\OneDrive - st.shou.edu.cn\\桌面\\Rcistarget_database\\target_gene_list.txt",sep="\t")
# motif Ranking files from https://resources.aertslab.org/cistarget/: 
motifRankings <- importRankings(dbFile = "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")

file_prefix <- "sc_cell_E14_P33_107kmeans_2cluster_"  # for output filenames
###############################################################################################################
geneList1 <- c("1110038B12Rik", "2410006H16Rik", "4632427E13Rik", "Aard", "Actl6a", "Aldh1a2", "Aldh2", "Anxa1", "Aplp2", "Apoc1", "Apod", "Apoe", "Atf6b", "Bax", "BC006965", "Bex3", "Bgn", "Blvrb", "Camp", "Camta1", "Casp3", "Cav2", "Ccl11", "Ccnd2", "Cct6a", "Cdc42ep3", "Cdk4", "Ciao2a", "Cldn11", "Cldn6", "Coch", "Col8a1", "Colec12", "Crym", "Ctsh", "Dazap1", "Ddr2", "Emx2", "Enpp2", "Ezr", "Fam43a", "Fgf20", "Fxyd1", "Gas5", "Gm14964", "Grb10", "Grn", "Gsn", "Gstt1", "H2-K1", "H2-T23", "H2afy", "Hes6", "Hmga2", "Hmgn1", "Hspb1", "Ier3", "Ifi30", "Igf2bp2", "Igfbp2", "Isg15", "Isl2", "Islr", "Itih2", "Lifr", "Lor", "Mdfi", "Meg3", "Mex3a", "Minar2", "Mrpl53", "Nbdy", "Nbl1", "Ndrg1", "Ntf3", "Nudt21", "Oat", "Pkm", "Prxl2a", "Psip1", "Ptn", "Rabac1", "Rarres2", "S100b", "Selenoh", "Sinhcaf", "Slc37a3", "Slco1c1", "Snu13", "Srsf6", "Sumo1", "Tac1", "Tapbp", "Tbx1", "Tcim", "Tead2", "Tent5a", "Tnfrsf11b", "Traf4", "Tuba1b", "Uchl3", "Vps72", "Vta1", "Wnt6", "Zfas1", "Zfp422", "Znrd2")
geneLists <- list(geneListName=geneList1)
# ---- Load gene sets to analyze (one set):
#geneset1 <- row.names(d) # get a vector of a gene \name
#geneLists <- list(sc107=geneset1)


# ---- Select motif database to use (i.e. organism and distance around TSS)
data(motifAnnotations_mgi) ## mouse: "motifAnnotations_mgi" ; human: "motifAnnotations_hgnc"
motifAnnotations_mgi[1:3]


# ---- key step:  Motif enrichment analysis:
motifEnrichmentTable_wGenes <- cisTarget(geneLists, 
                                         motifRankings = motifRankings,
                                         motifAnnot=motifAnnotations_mgi)
motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)  # adding logos


# ---- Output enriched motif results:
write.table(motifEnrichmentTable_wGenes,file = paste("6.RcisTarget_",file_prefix,".txt",sep=""), sep="\t",quote = F,row.names = F)


# ---- table visualization  (To get html file, click "Show in new window" in Viewer)
library(DT)  # BiocManager::install(c("DT"))
datatable(motifEnrichmentTable_wGenes_wLogo[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))



# ---- Build & Visualize candiate motif & target in network

##  get  candidate motif: 
candMotif <- motifEnrichmentTable_wGenes$motif[1:6]
candMotif <- "taipale_cyt_meth__KLF15_NCCMCGCCCMYN_FL"

## get connections for motif: 
incidenceMatrix <- getSignificantGenes(geneLists$geneListName, 
                                       motifRankings,
                                       signifRankingNames=candMotif,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix


## Generate edge and node files for network:
library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")

motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))


## Visualization:
library(visNetwork)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE)

## optional: exporting net file (could be open by Cytoscope)
write.table(edges,paste("6.RcisTarget_edges_",file_prefix,".txt",sep=""),quote = F,row.names = F)






