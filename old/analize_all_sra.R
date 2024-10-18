setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)

sra_geno<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")
recount3_metadata<-fread("/dcs04/hansen/data/recount_genotype/PCA/SRA/Recount3_metadata.tsv", header= T, sep = "\t",quote="")
recount3_metadata<-recount3_metadata[,c(2:5,163)]
sra_geno$seq_type<-recount3_metadata$seq_type[match(sra_geno$sample_id,recount3_metadata$external_id)]
sra_geno<-sra_geno[which(sra_geno$seq_type=="bulk"),]
SRA_metadata<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA_metadata.csv")
sra_geno$pred.type<-SRA_metadata$pred.type[match(sra_geno$sample_id,SRA_metadata$external_id)]
id<-which(sra_geno$pred.type=="rna-seq")

df<-data.frame(ref_ratio=double())
for(k in 1:nrow(sra_geno)){
  print(k)
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sra_geno$sample_id[k], "_ase.rda")))
  {
load(paste0("~/hansen_lab/ASE/test_ASE/", sra_geno$sample_id[k], "_ase.rda") ) #named ase_all or ase_df

if(exists("ase_all")){
      df[k,]<-median(ase_all$ref_ratio)
      rm(ase_all)
}else{
      df[k,]<-median(ase_df$ref_ratio)
      rm(ase_df)
    }
}}
#saveRDS(df,file= "~/plot/ASE/all_sra.rds")
df<-readRDS("~/plot/ASE/all_sra.rds")
df_bulk<-data.frame(ref_ratio=df$ref_ratio[id])
dim(df)
summary(as.numeric(df$ref_ratio),na.rm=T)
which(is.na(df$ref_ratio))

pdf(file="~/plot/ASE/allSRA_ref_ratio.pdf", width = 10, height = 4)
ggplot(df_bulk, aes(ref_ratio))+
  geom_histogram(bins = 100)+
  geom_vline(xintercept = 0.5, color="red")+
  geom_vline(xintercept = 0.55, color="purple")+
  geom_vline(xintercept = 0.45, color="blue")+
  labs(title="Density plot of ref-ratio for all the samples in SRA. [blue=0.45, red=0.5, purple=0.5]",
       subtitle= paste0("Only predicted bulk samples used: # of samples=", nrow(df_bulk)))

ggplot(df_bulk, aes(ref_ratio))+
  geom_density()+
  geom_vline(xintercept = 0.5, color="red")+
  geom_vline(xintercept = 0.55, color="purple")+
  geom_vline(xintercept = 0.45, color="blue")+
  labs(title="Density plot of ref-ratio for all the samples in SRA. [blue=0.45, red=0.5, purple=0.5]",
       subtitle= paste0("Only predicted bulk samples used: # of samples=", nrow(df_bulk)))
dev.off()

which(df[id,]>0.7)[3]

sra_geno[id,][which(sra_geno$study[id]=="ERP003467")[8],]
df[id,][1418]


df[id,][which(sra_geno$study[id]=="ERP003467")[8]]
