
setwd("~/ASE/")
library(tidyverse)
library(ggplot2)
library(ggridges)

battle<-read_tsv("data/battle.tsv")

sra_subset<-readRDS("~/test/m.rds")
#sra_subset<-readRDS("~/plot/ASE/sra_subset.rds")
#sra_met<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/metadata/all_SRA.csv")
sra_subset$cancer<-battle$cancer[match(sra_subset$study,battle$study)]
sra_subset$disease.category<-battle$disease.category[match(sra_subset$study,battle$study)]


table(sra_subset_nonCancer$cancer)
sra_subset_nonCancer<-sra_subset[-which(sra_subset$cancer=="cancer"),]
sra_subset_Cancer<-sra_subset[which(sra_subset$cancer=="cancer"),]

lib_good<-sra_subset_nonCancer[sra_subset_nonCancer$library_selection%in% c("cDNA","RT-PCR"),]
lib_bad<-sra_subset_nonCancer[-which(sra_subset_nonCancer$library_selection%in% c("cDNA","RT-PCR")),]
table(lib_bad$library_selection)
lib_bad<-lib_bad[order(lib_bad$library_selection),]
pdf(file="~/plot/ASE/library_selection_cDNA_RTPCT.pdf", width = 10, height = 6) #library_selection_NOTcDNA_RTPCT.pdf

data<-lib_good
for(k in 1:100){
  sample_id<-data$sample_id[k]
  study<-data$study[k]
  cc<-data$cancer[k]
  dd<-data$disease.category[k]
  ll<-data$library_selection[k]
  print(k)
  
  
  load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
  if(exists("ase_all")){
    ase_df<-ase_all
    rm(ase_all)
  } else{
  data$num_snp[k]<- nrow(ase_df)
  
  
  p0=ggplot(ase_df, aes(y=log10(alt), x=log10(ref)))+
    geom_point(alpha=0.4)+
    geom_smooth(method = lm, se = FALSE, color="blue")+
    geom_abline(intercept = 0, slope = 1,color="red")+
    geom_abline(intercept = 0.25, slope = 1,color="salmon")+
    geom_abline(intercept = -0.25, slope = 1,color="salmon")+
    labs(title= paste0(sample_id,"-", data$library_layout[k]),
         subtitle= paste0(study,"-",median(ase_df$ref_ratio),
                          "-", nrow(ase_df),"-",cc,"-", dd,"-", ll))
  print(p0)
  
  }
  
  rm(ase_df)
  
}
dev.off()
#-----------------------------------


good<-sra_subset_nonCancer %>% filter(ref_ratio==0.5)
bad<-sra_subset_nonCancer %>% filter(ref_ratio>=0.55)
low<-sra_subset_nonCancer %>% filter(ref_ratio<=0.48)

paired_noOverlap<-sra_subset_nonCancer %>% filter(library_layout=="paired",overlap>0)
paired_overlap<-sra_subset_nonCancer %>% filter(library_layout=="paired",overlap< (-26))

single_seqLong<-sra_subset_nonCancer %>% filter(library_layout=="single",seq_len>76)
single_seqShort<-sra_subset_nonCancer %>% filter(library_layout=="single",seq_len<50)

quantile(single_seqShort$num_snp,na.rm=T)

test_data<-sra_subset %>% filter(library_layout=="paired",overlap<17 & overlap>(-26))

test_data$num_snp<-NA
pdf(file="~/plot/ASE/single_seqShort.pdf", width = 10, height = 6)
data<-single_seqShort
for(k in 1:100){
  sample_id<-data$sample_id[k]
  study<-data$study[k]
  cc<-data$cancer[k]
  dd<-data$disease.category[k]
  base_mil<-data$num_bases_mil[k]
  print(k)
  
  
  load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
  if(exists("ase_all")){
    ase_df<-ase_all
    rm(ase_all)
  } else{
  data$num_snp[k]<- nrow(ase_df)
  
  
  p0=ggplot(ase_df, aes(y=log10(alt), x=log10(ref)))+
    geom_point(alpha=0.4)+
    geom_smooth(method = lm, se = FALSE, color="blue")+
    geom_abline(intercept = 0, slope = 1,color="red")+
    geom_abline(intercept = 0.25, slope = 1,color="salmon")+
    geom_abline(intercept = -0.25, slope = 1,color="salmon")+
    labs(title= paste0(sample_id,"-", data$library_layout[k]),
         subtitle= paste0(study,"-",round(median(ase_df$ref_ratio),3),
                          "-", nrow(ase_df),"-",cc,"-", dd, "-", base_mil))
  print(p0)
  
  
  rm(ase_df)
  
}}
dev.off()


