getwd()
install.packages('tidyverse')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install("enrichplot")
  library(ChIPseeker)
  library(clusterProfiler)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(DOSE)
  library(ReactomePA)
  library(ggplot2)
  library(DiffBind)
library(dplyr)
library(RMariaDB)
library(plyr)
library(readr)
library(enrichplot)
#getting my file frome my directory
#"E:/Admission/31724112/Indiana University/My Lab/ATAC-seq"
setwd("E:/Admission/31724112/Indiana University/My Lab/ATAC-seq/Peaks")
txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene
myfiles1<-list.files(pattern = "*.narrowpeak")
peakAnnoList<-lapply(myfiles1, annotatePeak, TxDb=txdb,
                 tssRegion=c(-3000, 3000),annoDb="org.Hs.eg.db")


#Making loop for peak annotations
myfiles2<-list(myfiles1[[1]],myfiles1[[2]],myfiles1[[3]],myfiles1[[4]],
               myfiles1[[5]])
myfiles3<-c("16hr_CIS", "CIS_RESIS", "DAC_CIS", "DAC_ONLY", 
              "UT")


for (i in (1:length(myfiles2))){
  n<-myfiles2[[i]]
  n<- annotatePeak(n, tssRegion=c(-3000, 3000),
                   TxDb=txdb, annoDb="org.Hs.eg.db")
  assign(paste(myfiles3[i],'Anno', sep=""), n)
}

#converting the Granges file to dataframe format
hr16_CIS<-as.data.frame(`16hr_CISAnno`)
CIS_RESIS<-as.data.frame(CIS_RESISAnno)
DAC_CIS<-as.data.frame(DAC_CISAnno)
DAC_ONLY<-as.data.frame(DAC_ONLYAnno)
UT<-as.data.frame(UTAnno)
#writting a table for the annotated genes
write.table(hr16_CISAnno, file ="hr16_CISAnno",sep='\t',
            row.names = FALSE)
write.table(CIS_RESISAnno, file ="CIS_RESISAnno",sep='\t',
            row.names = FALSE)
write.table(DAC_CISAnno, file ="DAC_CISAnno",sep='\t',
            row.names = FALSE)
write.table(DAC_ONLYAnno, file ="DAC_ONLYAnno",sep='\t',
            row.names = FALSE)
write.table(UTAnno, file ="UTAnno",sep='\t',
            row.names = FALSE)

#Visualize Genomic Annotation
pdf("hr16_CIS.pdf")
vennpie(`16hr_CISAnno`)
dev.off()

# doing peak files for enrichpathway
myfiles1
hr16_CIS_peak<-readPeakFile(myfiles1[[1]],as="GRanges")
CIS_RESIS_peak<-readPeakFile(myfiles1[[2]],as="GRanges")
DAC_CIS_peak<-readPeakFile(myfiles1[[3]],as="GRanges")
DAC_ONLY_peak<-readPeakFile(myfiles1[[4]],as="GRanges")
UT_peak<-readPeakFile(myfiles1[[5]],as="GRanges")
#Please note that you have to use data frame and peak file for
#functional enrichment analysis
#loop for enrichpathway in reactome package


mypeaks<-list(hr16_CIS_peak,CIS_RESIS_peak,DAC_CIS_peak,DAC_ONLY_peak,
              UT_peak)
mypeaks1<-c("hr16_CIS", "CIS_RESIS", "DAC_CIS", "DAC_ONLY", 
            "UT")
hr16_CIS.gene<-seq2gene(hr16_CIS_peak,tssRegion = c(-1000, 1000), 
                        flankDistance = 3000, TxDb = txdb)

for (i in (1:length(mypeaks))){
  n<-mypeaks[[i]]
  n<- seq2gene(n, tssRegion = c(-1000,1000), flankDistance = 3000, TxDb=txdb,
               sameStrand = FALSE)
  assign(paste(mypeaks1[i],'gene', sep="."), n)
}

path_16hrcis <- enrichPathway(hr16_CIS.gene, pvalueCutoff = 0.01)
path_CISresis<-enrichPathway(CIS_RESIS.gene, pvalueCutoff = 0.01)
path_DAC_CIS<-enrichPathway(DAC_CIS.gene, pvalueCutoff = 0.01)
path_DAC_only<-enrichPathway(DAC_ONLY.gene, pvalueCutoff = 0.01)
path_UT<-enrichPathway(UT.gene, pvalueCutoff = 0.01)

pdf("hr16_cis_pathway.pdf",width=17, height=9)
dotplot(path_16hrcis, showCategory=15, orderBy = "x")
dev.off()

pdf("cis_resis_pathway.pdf",width=17, height=9)
dotplot(path_CISresis, showCategory=15, orderBy = "x")
dev.off()
pdf("DAC_cis_pathway.pdf",width=17, height=9)
dotplot(path_DAC_CIS, showCategory=15, orderBy = "x")
dev.off()
pdf("DAC_only_pathway.pdf",width=17, height=9)
dotplot(path_DAC_only, showCategory=15, orderBy = "x")
dev.off()
pdf("UT_pathway.pdf",width=17, height=9)
dotplot(path_UT, showCategory=15, orderBy = "x")
dev.off()

#ChIP peak data set comparison
#Functional profiles comparison

#Profile of several ChIP peak data binding to TSS region
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
gc()
tagMatrixList <- lapply(myfiles1, getTagMatrix, windows=promoter)
pdf("Average_peaks_TSS.pdf",width=17, height=9)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,
            resample=500, facet="row")
dev.off()

#ATAC peak annotation comparision
myfiles1

peakAnnoList1<-peakAnnoList
names(peakAnnoList1)<-c("hr16_cisplatin","Cis_resistance",
                        "DAC_cisplatin","DAC_only", "UT")
pdf("genomic_annotation.pdf",width=17, height=9)
plotAnnoBar(peakAnnoList1)
dev.off()
# distance relative to TSS
pdf("distance_to_TSS.pdf",width=17, height=9)
plotDistToTSS(peakAnnoList1)
dev.off()

#Functional profiles comparison
genes<-lapply(peakAnnoList1, function(i) as.data.frame(i)$geneId)
names(genes)<-sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.01,
                           pAdjustMethod = "BH")
pdf("KEGG_pathway.pdf",width=20,height=9)
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
dev.off()

compGO <- compareCluster(geneCluster   = genes,
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.01,OrgDb='org.Hs.eg.db',ont='BP',
                           pAdjustMethod = "BH")
pdf("GO_BF_pathway.pdf",width=20,height=9)
dotplot(compGO, showCategory = 15, title = "Gene Ontology (BF) Analysis")
dev.off()

#Overlap of peaks and annotated gene
genes<-lapply(peakAnnoList1, function(i) as.data.frame(i)$geneId)

pdf("overlapped_peaks.pdf",width=20,height=9)
vennplot(genes)
dev.off()