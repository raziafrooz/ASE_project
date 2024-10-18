setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)
library(ggplot2)
library(ggridges)


sra_geno<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")
recount3_metadata<-fread("/dcs04/hansen/data/recount_genotype/PCA/SRA/Recount3_metadata.tsv", header= T, sep = "\t",quote="")
recount3_metadata<-recount3_metadata[,c(2:5,163)]
sra_geno$seq_type<-recount3_metadata$seq_type[match(sra_geno$sample_id,recount3_metadata$external_id)]
sra_geno<-sra_geno[which(sra_geno$seq_type=="bulk"),]
SRA_metadata<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA_metadata.csv")
sra_geno$pred.type<-SRA_metadata$pred.type[match(sra_geno$sample_id,SRA_metadata$external_id)]
id<-which(sra_geno$pred.type=="rna-seq")

df<-readRDS("~/plot/ASE/all_sra.rds")
df_bulk<-data.frame(ref_ratio=df$ref_ratio[id])
sra_geno_bulk<-sra_geno[id,]
df_bulk<-cbind(df_bulk,sra_geno_bulk[,2:3])


df_bulk$geno_cov<-NA
df_bulk$het_cov<-NA
df_bulk$homo_ref_cov<-NA
df_bulk$homo_alt_cov<-NA

for(k in 1:nrow(sra_geno_bulk)){
  print(k)
  sample_id<-sra_geno_bulk$sample_id[k]
  study<-sra_geno_bulk$study[k]
  geno_dir<-sra_geno_bulk$genotypedSamples[k]
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda")))
  {
    
    ase_df<-as_tibble(read.csv(geno_dir)) %>% 
      filter(coverage>=8)
    df_bulk$geno_cov[k]<-sum(ase_df$coverage)
    df_bulk$het_cov[k]<-sum(ase_df$coverage[ase_df$pred_genotype==2])
    df_bulk$homo_ref_cov[k]<-sum(ase_df$coverage[ase_df$pred_genotype==1])
    df_bulk$homo_alt_cov[k]<-sum(ase_df$coverage[ase_df$pred_genotype==3])
    
  }
}

df_bulk$dis<-"low"
df_bulk$dis[df_bulk$ref_ratio>0.52]<-"high"
df_bulk$dis[df_bulk$ref_ratio>=0.48 & df_bulk$ref_ratio<=0.52]<-"high"


pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
ggplot(plot_df2,aes(y=interval, x=nSNP, fill=avg_len_interv))+
  geom_density_ridges(alpha=0.5)+
  labs(title="Number of heterozygous SNPs")

