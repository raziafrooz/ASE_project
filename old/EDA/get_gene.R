setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)
library(ggplot2)
library(ggridges)

df<-readRDS("~/plot/ASE/sra_qc.rds")
df_nor<-readRDS("~/plot/ASE/sra_qc_normal.rds")
colnames(df)
plot1<-rbind(df[,c(1:3,5,6)], df_nor)



library(AnnotationHub)
ah <- AnnotationHub()
ah <- query(ah, c("v26","GENCODE","Homo sapiens","GRch38")) #v26 was used in GTExV8
TxDb<-ah[["AH75155"]]
#gene_grch38<-genes(TxDb,columns=c("TXID", "TXNAME"))
tx_38<-transcripts(TxDb,columns=c("TXNAME","GENEID"))



k=1
for(k in 102:nrow(plot1)){
  print(k)
  sample_id<-plot1$sample_id[k]
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda")))
  {
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda"))
    if(exists("ase_all")){
      ase_df<-ase_all
      rm(ase_all)
    }
    
    ase_gr<-makeGRangesFromDataFrame(ase_df, seqnames.field="chr",
                                     start.field="pos", end.field="pos")
    ov<-findOverlaps(tx_38,ase_gr)
    
    plot1$num_genes[k]<-length(unique(unlist(tx_38$GENEID[queryHits(ov)])))
  }
}
plot1[10:20,c(1,6)]
df_bulk$dis<-"low"
df_bulk$dis[df_bulk$ref_ratio>0.52]<-"high"
df_bulk$dis[df_bulk$ref_ratio>=0.48 & df_bulk$ref_ratio<=0.52]<-"high"


pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
ggplot(plot_df2,aes(y=interval, x=nSNP, fill=avg_len_interv))+
  geom_density_ridges(alpha=0.5)+
  labs(title="Number of heterozygous SNPs")

