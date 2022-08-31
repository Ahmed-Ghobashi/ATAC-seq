#loading SampleSheet data
#Diff. peak Binding is pweformed through conda environment 
getwd()

setwd("E:/Admission/31724112/Indiana University/My Lab/ATAC-seq/Peaks/Diff Peaks")
#get the files and convert datafram to Granges
myfil1<-list.files(pattern = "*.txt")
Cis_Resis_vs_DAC_Cis<-read.table("Cis_Resis_vs_DAC_Cis.txt",header=TRUE)
UT_vs_DAC_Cis_GR<-makeGRangesFromDataFrame(UT_vs_DAC_Cis,
                                                  keep.extra.columns=TRUE,
                                                  seqnames.field="seqnames",
                                                  start.field="start",
                                                  end.field="end")
Cis_Resis_vs_DAC<-read.table("Cis_Resis_vs_DAC.txt",header=TRUE) 
Cis_vs_Cis_resis<-read.table("Cis_vs_Cis_resis.txt",header=TRUE)
Cis_vs_DAC_CIS<-read.table("Cis_vs_DAC_CIS.txt",header=TRUE)
UT_vs_CIS<-read.table("UT_vs_CIS.txt",header=TRUE)
UT_vs_DAC<-read.table("UT_vs_DAC.txt",header=TRUE)
Cis_vs_DAC<-read.table("Cis_vs_DAC.txt",header=TRUE)
DAC_Cis_vs_DAC<-read.table("DAC_Cis_vs_DAC.txt",header=TRUE)
UT_vs_Cis_resis<-read.table("UT_vs_Cis_resis.txt",header=TRUE)
UT_vs_DAC_Cis<-read.table("UT_vs_DAC_Cis.txt",header=TRUE)

#Making loop for peak annotations
myfil2<-list(Cis_Resis_vs_DAC_Cis_GR,Cis_Resis_vs_DAC_GR,Cis_vs_Cis_resis_GR
             ,Cis_vs_DAC_CIS_GR,UT_vs_CIS_GR,UT_vs_DAC_GR,
         Cis_vs_DAC_GR,DAC_Cis_vs_DAC_GR,UT_vs_Cis_resis_GR,UT_vs_DAC_Cis_GR)
memory.limit(size=10000)             

myfil3<-c("Cis_Resis_vs_DAC_Cis", "Cis_Resis_vs_DAC",
          "Cis_vs_Cis_resis", "Cis_vs_DAC_CIS","UT_vs_CIS","UT_vs_DAC",
          "Cis_vs_DAC","DAC_Cis_vs_DAC","UT_vs_Cis_resis","UT_vs_DAC_Cis"
            )


for (i in (1:length(myfil2))){
  n<-myfil2[[i]]
  n<- annotatePeak(n, tssRegion=c(-3000, 3000),
                   TxDb=txdb, annoDb="org.Hs.eg.db")
  assign(paste(myfil3[i],'Annotat', sep="."), n)
}

#converting the Granges file to dataframe format
Cis_Resis_vs_DAC.Annotat1<-as.data.frame(Cis_Resis_vs_DAC.Annotat)
Cis_Resis_vs_DAC_Cis.Annotat1<-as.data.frame(Cis_Resis_vs_DAC_Cis.Annotat)
Cis_vs_Cis_resis.Annotat1<-as.data.frame(Cis_vs_Cis_resis.Annotat)
Cis_vs_DAC_CIS.Annotat1<-as.data.frame(Cis_vs_DAC_CIS.Annotat)
UT_vs_CIS.Annotat1<-as.data.frame(UT_vs_CIS.Annotat)
UT_vs_DAC.Annotat1<-as.data.frame(UT_vs_DAC.Annotat)
Cis_vs_DAC.Annotat1<-as.data.frame(Cis_vs_DAC.Annotat)
DAC_Cis_vs_DAC.Annotat1<-as.data.frame(DAC_Cis_vs_DAC.Annotat)
UT_vs_Cis_resis.Annotat1<-as.data.frame(UT_vs_Cis_resis.Annotat)
UT_vs_DAC_Cis.Annotat1<-as.data.frame(UT_vs_DAC_Cis.Annotat)
#writting a table for the annotated Diff. peaks
getwd()
setwd("E:/Admission/31724112/Indiana University/My Lab/ATAC-seq/Peaks/Annotated peaks/Diff. Peak list")

write.table(UT_vs_DAC_Cis.Annotat1, file ="UT_vs_DAC_Cis.Annotat1.txt",
            sep='\t',
            row.names = FALSE)

#next step is to perform functional analysis for UT_VS Cis, UT_VS_Cis_resis
#and Cis_vs_cis_resis
#select the significant value FDR< 0.05
UT_vs_CIS.Annotat1<-UT_vs_CIS.Annotat1%>% filter(FDR < 0.05)
UT_vs_Cis_resis.Annotat1<-UT_vs_Cis_resis.Annotat1%>% filter(FDR < 0.05)
Cis_vs_Cis_resis.Annotat1<-Cis_vs_Cis_resis.Annotat1%>% filter(FDR < 0.05)

UT_VS_Cis_UP<-UT_vs_CIS.Annotat1%>% filter(FDR < 0.05 & Fold > 0)
UT_VS_Cis_down<-UT_vs_CIS.Annotat1%>% filter(FDR < 0.05 & Fold < 0)
UT_VS_Cis_resis_up<-UT_vs_Cis_resis.Annotat1%>% filter(FDR < 0.05 & Fold > 0)
UT_VS_Cis_resis_down<-UT_vs_Cis_resis.Annotat1%>% filter(FDR < 0.05 & Fold < 0)
Cis_vs_Cis_resis_up<-Cis_vs_Cis_resis.Annotat1%>% filter(FDR < 0.05 & Fold > 0)
Cis_vs_Cis_resis_down<-Cis_vs_Cis_resis.Annotat1%>% filter(FDR < 0.05 & Fold < 0)
#write bed files to use for data visualization 

write.table(Cis_vs_Cis_resis_down, file="Cis_vs_Cis_resis_down.bed", sep="\t", quote=F,
            row.names=F, col.names=F)

##DATA visualization and functional analysis

#Prepare the input , gene list around TSS
UT_vs_CIS.list<-subset(UT_vs_CIS.Annotat1.TSS,select=c(geneId,Fold))
UT_vs_CIS.resis.list<-subset(UT_vs_Cis_resis.Annotat1.TSS,select=c(geneId,Fold))
Cis_vs_CIS.resis.list<-subset(Cis_vs_Cis_resis.Annotat1.TSS,select=c(geneId,Fold))

## assume that 1st column is ID
## 2nd column is fold change
UTCIS.List<-UT_vs_CIS.list[,2]
UTCISresis.List<-UT_vs_CIS.resis.list[,2]
CIS.CISresis.List<-Cis_vs_CIS.resis.list[,2]

names(UTCIS.List)<-as.character(UT_vs_CIS.list[,1])
names(UTCISresis.List)<-as.character(UT_vs_CIS.resis.list[,1])
names(CIS.CISresis.List)<-as.character(Cis_vs_CIS.resis.list[,1])

UTCIS.List<-sort(UTCIS.List,decreasing = TRUE)
UTCISresis.List<-sort(UTCISresis.List,decreasing = TRUE)
CIS.CISresis.List<-sort(CIS.CISresis.List,decreasing = TRUE)

#doing GSEA
UTCIS.GO.BP<-gseGO(geneList     = UTCIS.List,
      OrgDb        = org.Hs.eg.db,
      ont          = "BP",
      nPerm        = 1000,
      minGSSize    = 3,
      maxGSSize    = 800,
      pvalueCutoff = 0.05,
      verbose      = FALSE)#it require alot of memory

#doing KEGG analysis
kegg_organism = "hsa"
UTCIS.KEGG <- enrichKEGG(gene= names(UTCIS.List),
               organism     = kegg_organism,
               pvalueCutoff = 0.05)
UTCIS.KEGG<-enrichDGN(names(UTCIS.List))
UTCIS.GSEA<-gseNCG(UTCIS.List, nPerm=10000)

pdf("UT vs CIS.GSEA.pdf",width=15,height=9)
ridgeplot(UTCIS.GSEA)
dev.off()

#combine all Dataframe into one big dataframe

UT_CIS_CISresis<-list(UT_VS_Cis_UP.TSS,UT_VS_Cis_down.TSS,UT_VS_Cis_resis_up.TSS,
                         UT_VS_Cis_resis_down.TSS, Cis_vs_Cis_resis_up.TSS,
                         Cis_vs_Cis_resis_down.TSS)
names(UT_CIS_CISresis)<-c("UTvsCIS_UP","UTvsCIS_down","UTvsCISresis_up",
                        "UTvsCISresis_down","CISvsCISresis_up", 
                        "CISvsCISres_down")
UT_CIS_CISresis$UTvsCISresis_down

UT_CIS_CISresis.ALL<-list(UT_vs_CIS.Annotat1.TSS,UT_vs_Cis_resis.Annotat1.TSS,
                          Cis_vs_Cis_resis.Annotat1.TSS)

names(UT_CIS_CISresis.ALL)<-c("UTvsCIS",
                        "UTvsCIS.resis","CISvsCIS.resis" )

#functional profile analysis


genes2<-lapply(UT_CIS_CISresis.ALL, function(i) as.data.frame(i)$geneId)
genes1<-lapply(UT_CIS_CISresis, function(i) as.data.frame(i)$geneId)

pdf("ven.DIff.pdf")
vennplot(genes2)
dev.off()
names(genes1)<-sub("_", "\n", names(genes))
compKEGG2 <- compareCluster(geneCluster   = genes1,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
pdf("KEGG_pathway2.groups.TSS.pdf",width=15,height=9)
dotplot(compKEGG2, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
dev.off()

compGO2 <- compareCluster(geneCluster   = genes1,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,OrgDb='org.Hs.eg.db',ont='BP',
                         pAdjustMethod = "BH")
pdf("GO_BP_pathway.Groups.TSS.pdf",width=20,height=9)
dotplot(compGO2, showCategory = 15, title = "Gene Ontology (BP) Analysis")
dev.off()
#
###convert the datafram to GRanges
UT_VS_Cis_UP_peak<-readPeakFile(UT_VS_Cis_UP_peak,as="GRanges")
UT_VS_Cis_down_peak<-readPeakFile(UT_VS_Cis_down,as="GRanges")
UT_VS_Cis_resis_up_peak<-readPeakFile(UT_VS_Cis_resis_up,as="GRanges")
UT_VS_Cis_resis_down_peak<-readPeakFile(UT_VS_Cis_resis_down,as="GRanges")
Cis_vs_Cis_resis_up_peak<-readPeakFile(Cis_vs_Cis_resis_up,as="GRanges")
Cis_vs_Cis_resis_down_peak<-readPeakFile(Cis_vs_Cis_resis_down,as="GRanges")

Cis_vs_Cis_resis_down_peak<-makeGRangesFromDataFrame(Cis_vs_Cis_resis_down,
                         keep.extra.columns=TRUE,
                         seqnames.field="seqnames",
                         start.field="start",
                         end.field="end")


##making loop for seq2gene
#loop for enrichpathway in reactome package

mypeak<-list(UT_VS_Cis_UP_peak,UT_VS_Cis_down_peak,UT_VS_Cis_resis_up_peak
             ,UT_VS_Cis_resis_down_peak,Cis_vs_Cis_resis_up_peak,
             Cis_vs_Cis_resis_down_peak)
              
mypeak1<-c("UT_CIS.up","UT_CIS.down", "UT_CIS.RESIS.up",
           "UT_CIS.RESIS.down", "CIS_CIS.RESIS.up", "CIS_CIS.RESIS.down") 

for (i in (1:length(mypeak))){
  n<-mypeak[[i]]
  n<- seq2gene(n, tssRegion = c(-1000,1000), flankDistance = 3000, TxDb=txdb,
               sameStrand = FALSE)
  assign(paste(mypeak1[i],'gene', sep="."), n)
}
      
path_UTCIS_UP<-enrichPathway(UT_CIS.up.gene, pvalueCutoff = 0.05)
path_UTCIS_down<-enrichPathway(UT_CIS.down.gene, pvalueCutoff = 0.05)
path_UTRESIS_up<-enrichPathway(UT_CIS.RESIS.up.gene, pvalueCutoff = 0.05)
path_UTRESIS_down<-enrichPathway(UT_CIS.RESIS.down.gene, pvalueCutoff = 0.05)
path_CIS.CISRESIS_up<-enrichPathway(CIS_CIS.RESIS.up.gene,pvalueCutoff = 0.05)
path_CIS.CISRESIS_down<-enrichPathway(CIS_CIS.RESIS.down.gene,pvalueCutoff = 0.05)

pdf("path_UTCIS_UP.pdf",width=17, height=9)
dotplot(path_UTCIS_UP, showCategory=15, orderBy = "x",title="UT vs CIS_UP")
dev.off()

pdf("path_UTCIS_down.pdf",width=17, height=9)
dotplot(path_UTCIS_down, showCategory=15, orderBy = "x",title="UT vs CIS_down")
dev.off()
pdf("path_UTRESIS_up.pdf",width=17, height=9)
dotplot(path_UTRESIS_up, showCategory=15, orderBy = "x",
        title="UT vs cis.RESIS_up")
dev.off()
pdf("path_UTRESIS_down.pdf",width=17, height=9)
dotplot(path_UTRESIS_down, showCategory=15, orderBy = "x",
        title="UT vs Cis.RESIS_down")
dev.off()
pdf("path_CIS.CISRESIS_up.pdf",width=17, height=9)
dotplot(path_CIS.CISRESIS_up, showCategory=15, orderBy = "x",
        title="CIS vs CIS.RESIS_up")
dev.off()
pdf("path_CIS.CISRESIS_down.pdf",width=17, height=9)
dotplot(path_CIS.CISRESIS_down, showCategory=15, orderBy = "x",
        title="CIS vs CIS.RESIS_down")
dev.off()

#filter the genes around the TSS (+-3kb)

UT_vs_CIS.Annotat1.TSS<-filter(UT_vs_CIS.Annotat1,
                  UT_vs_CIS.Annotat1$annotation %in% c("Promoter (<=1kb)",
              "Promoter (1-2kb)" ,"Promoter (2-3kb)","Downstream (<1kb)",
                                "Downstream (1-2kb)","Downstream (2-3kb)"))
  
UT_vs_Cis_resis.Annotat1.TSS<-filter(UT_vs_Cis_resis.Annotat1,
             UT_vs_Cis_resis.Annotat1$annotation %in% c("Promoter (<=1kb)",
            "Promoter (1-2kb)" ,"Promoter (2-3kb)","Downstream (<1kb)",
                              "Downstream (1-2kb)","Downstream (2-3kb)"))
Cis_vs_Cis_resis.Annotat1.TSS<-filter(Cis_vs_Cis_resis.Annotat1,
           Cis_vs_Cis_resis.Annotat1$annotation %in% c("Promoter (<=1kb)",
            "Promoter (1-2kb)" ,"Promoter (2-3kb)","Downstream (<1kb)",
                          "Downstream (1-2kb)","Downstream (2-3kb)"))

UT_VS_Cis_UP.TSS<-UT_vs_CIS.Annotat1.TSS%>% filter(FDR < 0.05 & Fold > 0)
UT_VS_Cis_down.TSS<-UT_vs_CIS.Annotat1.TSS%>% filter(FDR < 0.05 & Fold < 0)
UT_VS_Cis_resis_up.TSS<-UT_vs_Cis_resis.Annotat1.TSS%>% filter(FDR < 0.05 & Fold > 0)
UT_VS_Cis_resis_down.TSS<-UT_vs_Cis_resis.Annotat1.TSS%>% filter(FDR < 0.05 & Fold < 0)
Cis_vs_Cis_resis_up.TSS<-Cis_vs_Cis_resis.Annotat1.TSS%>% filter(FDR < 0.05 & Fold > 0)
Cis_vs_Cis_resis_down.TSS<-Cis_vs_Cis_resis.Annotat1.TSS%>% filter(FDR < 0.05 & Fold < 0)

write.table(Cis_vs_Cis_resis_down.TSS, file ="Cis_vs_Cis_resis_down.TSS.txt",
            sep='\t',
            row.names = FALSE)

##DATA visualization and functional analysis

#Prepare the input
UT_VS_Cis_UP.TSS.list<-subset(UT_VS_Cis_UP.TSS,select=c(geneId,Fold))
UT_VS_Cis_down.TSS.list<-subset(UT_VS_Cis_down.TSS,select=c(geneId,Fold))
UT_VS_Cis_resis_up.TSS.list<-subset(UT_VS_Cis_resis_up.TSS,select=c(geneId,Fold))
UT_VS_Cis_resis_down.TSS.list<-subset(UT_VS_Cis_resis_down.TSS,select=c(geneId,Fold))
Cis_vs_Cis_resis_up.TSS.list<-subset(Cis_vs_Cis_resis_up.TSS,select=c(geneId,Fold))
Cis_vs_Cis_resis_down.TSS.list<-subset(Cis_vs_Cis_resis_down.TSS,select=c(geneId,Fold))

## assume that 1st column is ID
## 2nd column is fold change
UTCIS.up.List<-UT_VS_Cis_UP.TSS.list[,2]
UTCIS.down.List<-UT_VS_Cis_down.TSS.list[,2]
UTCISresis.up.List<-UT_VS_Cis_resis_up.TSS.list[,2]
UTCISresis.down.List<-UT_VS_Cis_resis_down.TSS.list[,2]
CIS.CISresis.up.List<-Cis_vs_Cis_resis_up.TSS.list[,2]
CIS.CISresis.down.List<-Cis_vs_Cis_resis_down.TSS.list[,2]

names(UTCIS.up.List)<-as.character(UT_VS_Cis_UP.TSS.list[,1])
names(UTCIS.down.List)<-as.character(UT_VS_Cis_down.TSS.list[,1])
names(UTCISresis.up.List)<-as.character(UT_VS_Cis_resis_up.TSS.list[,1])
names(UTCISresis.down.List)<-as.character(UT_VS_Cis_resis_down.TSS.list[,1])
names(CIS.CISresis.up.List)<-as.character(Cis_vs_Cis_resis_up.TSS.list[,1])
names(CIS.CISresis.down.List)<-as.character(Cis_vs_Cis_resis_down.TSS.list[,1])

UTCIS.up.List<-sort(UTCIS.up.List,decreasing = TRUE)
UTCIS.down.List<-sort(UTCIS.down.List,decreasing = TRUE)
UTCISresis.up.List<-sort(UTCISresis.up.List,decreasing = TRUE)
UTCISresis.down.List<-sort(UTCISresis.down.List,decreasing = TRUE)
CIS.CISresis.up.List<-sort(CIS.CISresis.up.List,decreasing = TRUE)
CIS.CISresis.down.List<-sort(CIS.CISresis.down.List,decreasing = TRUE)

#using H for Gene Ontology(hallmark) like apoptosis, WNT (Hallmark analysis)

m_t4g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

em4.UTCIS.up <- enricher(names(UTCIS.up.List) , TERM2GENE=m_t4g)
em4.UTCIS.down <- enricher(names(UTCIS.down.List) , TERM2GENE=m_t4g)
em4.UTCISresis.up <- enricher(names(UTCISresis.up.List) , TERM2GENE=m_t4g)
em4.UTCISresis.down <- enricher(names(UTCISresis.down.List) , TERM2GENE=m_t4g)
em4.CIS.CISresis.up <- enricher(names(CIS.CISresis.up.List) , TERM2GENE=m_t4g)
em4.CIS.CISresis.down <- enricher(names(CIS.CISresis.down.List) , TERM2GENE=m_t4g)



pdf("H.UTCIS.up.barplot.pdf",height=7,width=10)
barplot(em4.UTCIS.up,showCategory=20, main="UT vs CIS.UP")
dev.off()
pdf("H.UTCIS.down.barplot.pdf",height=7,width=10)
barplot(em4.UTCIS.down,showCategory=20, main="UT vs CIS.down")
dev.off()
pdf("H.UTCISresis.up.barplot.pdf",height=7,width=10)
barplot(em4.UTCISresis.up,showCategory=20, main="UT vs CIS.Resis.up")
dev.off()
pdf("H.UTCISresis.down.barplot.pdf",,height=7,width=10)
barplot(em4.UTCISresis.down,showCategory=20, main="UT vs CIS.Resis.down")
dev.off()
pdf("H.CIS.CISresis.up.barplot.pdf",width=10, height=7)
barplot(em4.CIS.CISresis.up,showCategory=20, main="CIS vs CIS.Resis.up")
dev.off()

#using Biological process (BP)of C5 for Gene Ontology
m_t7g <- msigdbr(species = "Homo sapiens", category = "C5",subcategory="BP") %>% 
  dplyr::select(gs_name, entrez_gene)
dev.off()


em7.UTCIS.up <- enricher(names(UTCIS.up.List) , TERM2GENE=m_t7g)
em7.UTCIS.down <- enricher(names(UTCIS.down.List) , TERM2GENE=m_t7g)
em7.UTCISresis.up <- enricher(names(UTCISresis.up.List) , TERM2GENE=m_t7g)
em7.UTCISresis.down <- enricher(names(UTCISresis.down.List) , TERM2GENE=m_t7g)
em7.CIS.CISresis.up <- enricher(names(CIS.CISresis.up.List) , TERM2GENE=m_t7g)
em7.CIS.CISresis.down <- enricher(names(CIS.CISresis.down.List) , TERM2GENE=m_t7g)



pdf("BP.UTCIS.up.barplot.pdf",height=7,width=10)
barplot(em7.UTCIS.up,showCategory=20, main="UT vs CIS.UP")
dev.off()
pdf("BP.UTCIS.down.barplot.pdf",height=7,width=10)
barplot(em7.UTCIS.down,showCategory=20, main="UT vs CIS.down")
dev.off()
pdf("BP.UTCISresis.up.barplot.pdf",height=7,width=10)
barplot(em7.UTCISresis.up,showCategory=20, main="UT vs CIS.Resis.up")
dev.off()
pdf("BP.UTCISresis.down.barplot.pdf",,height=7,width=10)
barplot(em7.UTCISresis.down,showCategory=20, main="UT vs CIS.Resis.down")
dev.off()
pdf("BP.CIS.CISresis.up.barplot.pdf",width=10, height=7)
barplot(em7.CIS.CISresis.up,showCategory=20, main="CIS vs CIS.Resis.up")
dev.off()
pdf("BP.CIS.CISresis.down.barplot.pdf",width=10, height=7)
barplot(em7.CIS.CISresis.down,showCategory=20,main="CIS vs CIS.Resis.down" )
dev.off()

#using Biological process (BP)of C5 for Gene Ontology
m_t8g <- msigdbr(species = "Homo sapiens", category = "C5",subcategory="MF") %>% 
  dplyr::select(gs_name, entrez_gene)



em8.UTCIS.up <- enricher(names(UTCIS.up.List) , TERM2GENE=m_t8g)
em8.UTCIS.down <- enricher(names(UTCIS.down.List) , TERM2GENE=m_t8g)
em8.UTCISresis.up <- enricher(names(UTCISresis.up.List) , TERM2GENE=m_t8g)
em8.UTCISresis.down <- enricher(names(UTCISresis.down.List) , TERM2GENE=m_t8g)
em8.CIS.CISresis.up <- enricher(names(CIS.CISresis.up.List) , TERM2GENE=m_t8g)
em8.CIS.CISresis.down <- enricher(names(CIS.CISresis.down.List) , TERM2GENE=m_t8g)



pdf("MF.UTCIS.up.barplot.pdf",height=7,width=10)
barplot(em8.UTCIS.up,showCategory=20, main="UT vs CIS.UP")
dev.off()
pdf("MF.UTCIS.down.barplot.pdf",height=7,width=10)
barplot(em8.UTCIS.down,showCategory=20, main="UT vs CIS.down")
dev.off()
pdf("MF.UTCISresis.up.barplot.pdf",height=7,width=10)
barplot(em8.UTCISresis.up,showCategory=20, main="UT vs CIS.Resis.up")
dev.off()
pdf("MF.UTCISresis.down.barplot.pdf",,height=7,width=10)
barplot(em8.UTCISresis.down,showCategory=20, main="UT vs CIS.Resis.down")
dev.off()
pdf("MF.CIS.CISresis.up.barplot.pdf",width=10, height=7)
barplot(em8.CIS.CISresis.up,showCategory=20, main="CIS vs CIS.Resis.up")
dev.off()
pdf("MF.CIS.CISresis.down.barplot.pdf",width=10, height=7)
barplot(em8.CIS.CISresis.down,showCategory=20,main="CIS vs CIS.Resis.down" )
dev.off()

memory.limit(size = 12000)
#using Biological process (KEGG pathway)of C2 for Gene Ontology (pathway analysis)
m_t6g <- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:KEGG") %>% 
  dplyr::select(gs_name, entrez_gene)

em6.UTCIS.up <- enricher(names(UTCIS.up.List) , TERM2GENE=m_t6g)
em6.UTCIS.down <- enricher(names(UTCIS.down.List) , TERM2GENE=m_t6g)
em6.UTCISresis.up <- enricher(names(UTCISresis.up.List) , TERM2GENE=m_t6g)
em6.UTCISresis.down <- enricher(names(UTCISresis.down.List) , TERM2GENE=m_t6g)
em6.CIS.CISresis.up <- enricher(names(CIS.CISresis.up.List) , TERM2GENE=m_t6g)
em6.CIS.CISresis.down <- enricher(names(CIS.CISresis.down.List) , TERM2GENE=m_t6g)



pdf("KEGG.UTCIS.up.dotplot.pdf",height=7,width=10)
dotplot(em6.UTCIS.up,showCategory=20, title="UT vs CIS.UP")
dev.off()
pdf("KEGG.UTCIS.down.dotplot.pdf",height=7,width=10)
dotplot(em6.UTCIS.down,showCategory=20, title="UT vs CIS.down")
dev.off()
pdf("MF.UTCISresis.up.dotplot.pdf",height=7,width=10)
dotplot(em6.UTCISresis.up,showCategory=20, title="UT vs CIS.Resis.up")
dev.off()
pdf("MF.UTCISresis.down.dotplot.pdf",,height=7,width=10)
dotplot(em6.UTCISresis.down,showCategory=20, title="UT vs CIS.Resis.down")
dev.off()
pdf("MF.CIS.CISresis.up.dotplot.pdf",width=10, height=7)
dotplot(em6.CIS.CISresis.up,showCategory=20, title="CIS vs CIS.Resis.up")
dev.off()
pdf("KEGG.CIS.CISresis.down.dotplot.pdf",width=10, height=7)
dotplot(em6.CIS.CISresis.down,showCategory=20,title="CIS vs CIS.Resis.down" )


#using Biological process (REACTOME pathway)of C2 for Gene Ontology (pathway analysis)
m_t11g <- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:REACTOME") %>% 
  dplyr::select(gs_name, entrez_gene)


em11.UTCIS.up <- enricher(names(UTCIS.up.List) , TERM2GENE=m_t11g)
em11.UTCIS.down <- enricher(names(UTCIS.down.List) , TERM2GENE=m_t11g)
em11.UTCISresis.up <- enricher(names(UTCISresis.up.List) , TERM2GENE=m_t11g)
em11.UTCISresis.down <- enricher(names(UTCISresis.down.List) , TERM2GENE=m_t11g)
em11.CIS.CISresis.up <- enricher(names(CIS.CISresis.up.List) , TERM2GENE=m_t11g)
em11.CIS.CISresis.down <- enricher(names(CIS.CISresis.down.List) , TERM2GENE=m_t11g)



pdf("React.UTCIS.up.dotplot.pdf",height=15,width=50)
dotplot(em11.UTCIS.up,showCategory=20, title="UT vs CIS.UP")
dev.off()
pdf("React.UTCIS.down.dotplot.pdf",height=7,width=10)
dotplot(em11.UTCIS.down,showCategory=20, title="UT vs CIS.down")
dev.off()
pdf("React.UTCISresis.up.dotplot.pdf",height=7,width=10)
dotplot(em11.UTCISresis.up,showCategory=20, title="UT vs CIS.Resis.up")
dev.off()
pdf("React.UTCISresis.down.dotplot.pdf",,height=7,width=10)
dotplot(em11.UTCISresis.down,showCategory=20, title="UT vs CIS.Resis.down")
dev.off()
pdf("React.CIS.CISresis.up.dotplot.pdf",width=10, height=7)
dotplot(em11.CIS.CISresis.up,showCategory=20, title="CIS vs CIS.Resis.up")
dev.off()
pdf("React.CIS.CISresis.down.dotplot.pdf",width=10, height=7)
dotplot(em11.CIS.CISresis.down,showCategory=20,title="CIS vs CIS.Resis.down" )

#using Biological process (PID pathway)of C2 for Gene Ontology (pathway analysis)
m_t12g <- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:BIOCARTA") %>% 
  dplyr::select(gs_name, entrez_gene)


em12.UTCIS.up <- enricher(names(UTCIS.up.List) , TERM2GENE=m_t12g)
em12.UTCIS.down <- enricher(names(UTCIS.down.List) , TERM2GENE=m_t12g)
em12.UTCISresis.up <- enricher(names(UTCISresis.up.List) , TERM2GENE=m_t12g)
em12.UTCISresis.down <- enricher(names(UTCISresis.down.List) , TERM2GENE=m_t12g)
em12.CIS.CISresis.up <- enricher(names(CIS.CISresis.up.List) , TERM2GENE=m_t12g)
em12.CIS.CISresis.down <- enricher(names(CIS.CISresis.down.List) , TERM2GENE=m_t12g)



pdf("Biocarta.UTCIS.up.dotplot.pdf",height=7,width=10)
dotplot(em12.UTCIS.up,showCategory=20, title="UT vs CIS.UP")
dev.off()
pdf("Biocarta.UTCIS.down.dotplot.pdf",height=7,width=10)
dotplot(em12.UTCIS.down,showCategory=20, title="UT vs CIS.down")
dev.off()
pdf("biocarta.UTCISresis.up.dotplot.pdf",height=7,width=10)
dotplot(em12.UTCISresis.up,showCategory=20, title="UT vs CIS.Resis.up")
dev.off()
pdf("Biocarta.UTCISresis.down.dotplot.pdf",,height=7,width=10)
dotplot(em12.UTCISresis.down,showCategory=20, title="UT vs CIS.Resis.down")
dev.off()
pdf("Biocarta.CIS.CISresis.up.dotplot.pdf",width=10, height=7)
dotplot(em12.CIS.CISresis.up,showCategory=20, title="CIS vs CIS.Resis.up")
dev.off()
pdf("Biocarta.CIS.CISresis.down.dotplot.pdf",width=10, height=7)
dotplot(em12.CIS.CISresis.down,showCategory=20,title="CIS vs CIS.Resis.down" )

#doing GSEA (in anacoda)

edo2.UTCIS.up <- gseNCG(UTCIS.up.List, nPerm=10000)
edo2.UTCIS.down <- gseNCG(UTCIS.down.List, nPerm=10000)
edo2.UTCISresis.up <- gseNCG(UTCISresis.up.List, nPerm=10000)
edo2.UTCISresis.down <- gseNCG(UTCISresis.down.List, nPerm=10000)
edo2.CIS.CISresis.up <- gseNCG(CIS.CISresis.up.List, nPerm=10000)
edo2.CIS.CISresis.down <- gseNCG(CIS.CISresis.down.List, nPerm=10000)
