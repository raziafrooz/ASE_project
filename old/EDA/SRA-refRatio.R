setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(recount3)
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

id_low<-which(df_bulk<0.48)
low_dist<-df_bulk[id_low,]
id_high<-which(df_bulk>0.52 & df_bulk<0.6 )
high_dist<-df_bulk[id_high,]
length(id_low)
#First look to see if there are any differences between the read depth between the two distribution
rm(sra_geno_bulk_high)
plot_df<-low_dist
plot_df$dis<-"low"
#make plot df for high dist:
high_dist$dis<-"high"
plot_df<-rbind(plot_df,high_dist)
dim(plot_df)
plot_df$numberOfBases<-NA
plot_df$avg_len<-NA
study<-"DRR013762"
for(study in unique(plot_df$study)){
  print(study)
url<-locate_url(
  study,
  "data_sources/sra",
  type = "metadata")

xx <-utils::read.delim(file_retrieve(url[4], verbose = FALSE))
plot_df$numberOfBases[which(plot_df$study == study)]<- xx$X.bases[na.omit(match(plot_df$sample_id,xx$external_id))]
plot_df$avg_len[which(plot_df$study == study)]<- xx$avg_len[na.omit(match(plot_df$sample_id,xx$external_id))]

}


xx[1,3]
plot_df[which(plot_df$study == "ERP011097"),]
library(ggridges)


plot_df<-plot_df %>% group_by(dis) %>% mutate(interval= cut_number(ref_ratio, 6)) %>% ungroup()

pdf(file="~/plot/ASE/SRA_qc.pdf", width = 10, height = 4)
ggplot(plot_df2,aes(y=interval, x=numberOfBases, fill=dis))+
  geom_density_ridges()+
  labs(title="Number of seq base between two groups.",
       subtitle= "colored based on high (ref_ratio>0.52) and low (ref_ratio<0.48) distribution",
       y="ref ratio category")+
  xlim(c(0,2e+10))

ggplot(plot_df2,aes(y=interval, x=avg_len, fill=dis))+
  geom_density_ridges()+
  labs(title="Average sequence length between the two groups",
      subtitle= "colored based on high (ref_ratio>0.52) and low (ref_ratio<0.48) distribution",
      y="ref ratio category",
      x="average seq length")+
  xlim(c(20,200))

dev.off()

plot_df2 %>% filter(dis=="high")
plot_df2[1000,]
plot_df2$sample_id[1]
k=841
plot_df2$star_unique_mapped<-NA
for(k in 1637:length(unique(plot_df2$study))){
  study<-unique(plot_df2$study)[k]
  print(k)
  url<-locate_url(
    study,
    "data_sources/sra",
    type = "metadata")
  
  xx <-utils::read.delim(file_retrieve(url[3], verbose = FALSE))
  plot_df2$multi_map[which(plot_df2$study == study)]<- xx$star.._of_reads_mapped_to_multiple_loci[na.omit(match(plot_df2$sample_id,xx$external_id))]
  plot_df2$star_all_mapped[which(plot_df2$study == study)]<- xx$star.all_mapped_reads[na.omit(match(plot_df2$sample_id,xx$external_id))]
  plot_df2$star_unique_mapped[which(plot_df2$study == study)]<- xx$star.uniquely_mapped_reads_number[na.omit(match(plot_df2$sample_id,xx$external_id))]
  
}



for(k in 1662:nrow(plot1)){
  print(k)
  sample_id<-plot1$sample_id[k]
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda")))
  {
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda"))
    if(exists("ase_all")){
      ase_df<-ase_all
      rm(ase_all)
    }
    
    
    plot1$ase_sig[k]<-sum(ase_df$q_val<=0.01)
    plot1$nSNP[k]<-nrow(ase_df)
    rm(ase_df)
  }
}


plot_df2<-plot_df2 %>% group_by(dis) %>% mutate(avg_len_interv= cut_number(avg_len, 2)) %>% ungroup()
plot_df2$avg_len_interv
pdf(file="~/plot/ASE/SRA_qc3.pdf", width = 10, height = 4)
ggplot(plot_df2,aes(y=interval, x=multi_map, fill=avg_len_interv))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,50))+
  labs(title="STAR:Number of reads mapped to multiple loci divided by number of input reads",
       x="% mapped to multiple loci")

ggplot(plot_df2,aes(y=interval, x=star_all_mapped, fill=avg_len_interv))+
  geom_density_ridges(alpha=0.5)+
  labs(title="STAR:Total number of reads aligned")

ggplot(plot_df2,aes(y=interval, x=star_unique_mapped, fill=avg_len_interv))+
  geom_density_ridges(alpha=0.5)+
  labs(title="STAR:Number of reads which mapped to a single locus")
dev.off() 

plot_df2$library_layout<-NA
for(k in 1:length(unique(plot_df2$study))){
  study<-unique(plot_df2$study)[k]
  print(k)
  url<-locate_url(
    study,
    "data_sources/sra",
    type = "metadata")
  
  xx <-utils::read.delim(file_retrieve(url[1], verbose = FALSE))
  plot_df2$library_layout[which(plot_df2$study == study)]<- xx$library_layout[na.omit(match(plot_df2$sample_id,xx$external_id))]
  
}
pdf(file="~/plot/ASE/SRA_qc2.pdf", width = 10, height = 4)
ggplot(plot_df2,aes(y=interval, x=avg_len, fill=library_layout))+
  geom_density_ridges(alpha=0.5)+
  labs(title="sequence average length for ref-ratio intervals based on library layout")
dev.off() 

colnames(plot_df2)

pdf(file="~/plot/ASE/SRA_qc4.pdf", width = 10, height = 4)
ggplot(plot_df2,aes(y=interval, x=nSNP, fill=avg_len_interv))+
  geom_density_ridges(alpha=0.5)+
  labs(title="Number of heterozygous SNPs")

ggplot(plot_df2,aes(y=interval, x=total50p, fill=avg_len_interv))+
  geom_density_ridges(alpha=0.5)+
  labs(title="Total read count for 50% of the data")

ggplot(plot_df2,aes(y=interval, x=ref50p, fill=avg_len_interv))+
  geom_density_ridges(alpha=0.5)+
  labs(title="ref read count for 50% of the data")

dev.off() 

pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
ggplot(plot_df,aes(numberOfBases))+
  geom_density()
dev.off()
plot_df$numberOfBases


#saveRDS(plot_df2, file="~/plot/ASE/sra_qc.rds")
