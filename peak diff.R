getwd()
#go to the annotated peak file
setwd("Annotated peaks")
myfiles_annot<-list.files(pattern = "*.CSV", all.files = FALSE)
list(myfiles_annot[[2]])
CIS_RESIS_Anno<-read.table("CIS_RESISAnno", header=TRUE)
DAC_CIS_Anno<-read.table("DAC_CISAnno", header=TRUE)
DAC_only_Anno<-read.table("DAC_ONLYAnno",header=TRUE)
hr16_CIS_Anno<-read.csv("hr16_CISAnno.csv")
UT_Anno<-read.table("UTAnno",header = TRUE)

CIS_RESIS_Anno<-as.data.frame(CIS_RESIS_Anno)
DAC_CIS_Anno<-as.data.frame(DAC_CIS_Anno)
DAC_only_Anno<-as.data.frame(DAC_only_Anno)
hr16_CIS_Anno<-as.data.frame(hr16_CIS_Anno)
UT_Anno<-as.data.frame(UT_Anno)

       

#making loop for getting genes around the promoter

filter(CIS_RESIS_Anno, CIS_RESIS_Anno$annotation %in% c("Promoter (<=1kb)",
              "Promoter (1-2kb)" ,"Promoter (2-3kb)","Downstream (<1kb)",
     "Downstream (1-2kb)","Downstream (2-3kb)"))

myXfiles2<-list(CIS_RESIS_Anno,DAC_CIS_Anno,DAC_only_Anno,hr16_CIS_Anno,
               UT_Anno)
myXfiles3<-c("CIS_RESIS", "DAC_CIS", "DAC_only", "hr16_CIS", 
            "UT")


for (i in (1:length(myXfiles2))){
  n<-myXfiles2[[i]]
  n<- filter(n, n$annotation %in% c("Promoter (<=1kb)",
       "Promoter (1-2kb)" ,"Promoter (2-3kb)","Downstream (<1kb)",
               "Downstream (1-2kb)","Downstream (2-3kb)"))
  assign(paste(myXfiles3[i],'TSS', sep="."), n)
}
#Remove the duplicates
CIS_RESIS.TSS<-CIS_RESIS.TSS[!duplicated(CIS_RESIS.TSS$SYMBOL),]
DAC_CIS.TSS<-DAC_CIS.TSS[!duplicated(DAC_CIS.TSS$SYMBOL),]
DAC_only.TSS<-DAC_only.TSS[!duplicated(DAC_only.TSS$SYMBOL),]
hr16_CIS.TSS<-hr16_CIS.TSS[!duplicated(hr16_CIS.TSS$SYMBOL),]
UT.TSS<-UT.TSS[!duplicated(UT.TSS$SYMBOL),]
#remove NA value
UT.TSS_gene<-na.omit(UT.TSS$SYMBOL)
CIS_RESIS.TSS_gene<-na.omit(CIS_RESIS.TSS$SYMBOL)
DAC_CIS.TSS_gene<-na.omit(DAC_CIS.TSS$SYMBOL)
DAC_only.TSS_gene<-na.omit(DAC_only.TSS$SYMBOL)
hr16_CIS.TSS_gene<-na.omit(hr16_CIS.TSS$SYMBOL)
UT_specific_gene<-na.omit(UT_specific_genes$SYMBOL)


#the plan is to get the genes specific for UT condition
UT_hr16<-anti_join(UT.TSS,hr16_CIS.TSS, by="SYMBOL")
UT_hr16_DAC<-anti_join(UT_hr16,DAC_only.TSS, by="SYMBOL")
UT_HR_DAC_CIS<-anti_join(UT_hr16_DAC,DAC_CIS.TSS,by="SYMBOL")
UT_specific_genes<-anti_join(UT_HR_DAC_CIS,CIS_RESIS.TSS,by="SYMBOL")

write.table(UT_specific_genes,file="UT_specific",sep="\t",row.names = TRUE,
          col.names = TRUE )
write.table(UT_hr16,file="UT_no16hr",sep="\t",row.names = TRUE,
            col.names = TRUE )

#getting 16hr Cis specific genes
hr16_UT<-anti_join(hr16_CIS.TSS, UT.TSS,by="SYMBOL")
hr16_UT_DAC<-anti_join(hr16_UT,DAC_only.TSS, by="SYMBOL")
HR16_UT_DAC_CIS<-anti_join(hr16_UT_DAC,DAC_CIS.TSS,by="SYMBOL")
hr16_specific_genes<-anti_join(HR16_UT_DAC_CIS,CIS_RESIS.TSS,by="SYMBOL")

write.table(hr16_specific_genes,file="16HR_specific",sep="\t",row.names = TRUE,
            col.names = TRUE )

#getting Cis-resis specific genes

cis_resis_UT<-anti_join(CIS_RESIS.TSS, UT.TSS,by="SYMBOL")
resis_UT_DAC<-anti_join(cis_resis_UT,DAC_only.TSS, by="SYMBOL")
resis_UT_DAC_CIS<-anti_join(resis_UT_DAC,DAC_CIS.TSS,by="SYMBOL")
cis_resis_specific_genes<-anti_join(resis_UT_DAC_CIS,hr16_CIS.TSS,by="SYMBOL")

write.table(cis_resis_specific_genes,
            file="cis_Resis_specific",sep="\t",row.names = TRUE,
            col.names = TRUE )

#getting DAC-cis specific

cis_resis_UT1<-anti_join(DAC_CIS.TSS,CIS_RESIS.TSS,by="SYMBOL")
resis_UT_DAC1<-anti_join(cis_resis_UT1,UT.TSS, by="SYMBOL")
resis_UT_DAC_CIS1<-anti_join(resis_UT_DAC1,DAC_only.TSS,by="SYMBOL")
DAC_cis_specific_genes<-anti_join(resis_UT_DAC_CIS1,hr16_CIS.TSS,by="SYMBOL")

write.table(DAC_cis_specific_genes,
            file="DAC_cis_specific",sep="\t",row.names = TRUE,
            col.names = TRUE )

#getting DAC only

cis_resis_UT2<-anti_join(DAC_only.TSS,DAC_CIS.TSS,by="SYMBOL")
resis_UT_DAC2<-anti_join(cis_resis_UT2,UT.TSS, by="SYMBOL")
resis_UT_DAC_CIS2<-anti_join(resis_UT_DAC2,CIS_RESIS.TSS,by="SYMBOL")
DAC_only_specific_genes<-anti_join(resis_UT_DAC_CIS2,hr16_CIS.TSS,by="SYMBOL")

write.table(DAC_only_specific_genes,
            file="DAC_only_specific",sep="\t",row.names = TRUE,
            col.names = TRUE )
#Profile of several ChIP peak data binding to TSS region
str(UT_specific_genes)
TxDb_UT<-makeTxDbFromUCSC(genome="hg38",tablename="knownGene",
                          transcript_ids=UT_specific_genes$SYMBOL)
UT_specific_peaks<-readPeakFile(UT_specific_genes,as="GRanges")
promoter_UT <- getPromoters(TxDb= TxDb_UT,
                         upstream=3000, downstream=3000)
keytypes(txdb)
tagMatrixList_UT <- lapply(myfiles1, getTagMatrix, windows=promoter_UT)
pdf("Average_peaks_TSS.pdf",width=17, height=9)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,
            resample=500, facet="row")
dev.off()


library(VennDiagramn)
#making venndigram
venn.diagram(x=list(CIS_RESIS.TSS_gene,DAC_CIS.TSS_gene,
    hr16_CIS.TSS_gene,UT.TSS_gene,UT_specific_gene),
    category.names =c("CIS_RESIS","DAC_CIS","16hr_cis","UT",
      "UT_specific"),filename = "venndigram1.png",output=TRUE )

#diffBinding

myfiles1_ANNOTAT<-list.files(pattern="*.csv")
myfiles1_ANNOTAT1<-ldply(myfiles1_ANNOTAT,read_csv)
read.csv()
peakAnnoList1
diff_binding<-dba(sampleSheet =myfiles1_ANNOTAT1)

