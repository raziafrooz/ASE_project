setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(ggplot2)

study<-"Lung"

study_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rds"))

#--------------------------------------------
#get true gtex:
met<-read.csv("data/GTEx_metadata.csv")
tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

i<- which(tissues_names$full_name==study)
print(tissues_names$full_name[i])
abb<-tissues_names$abb[i]
tissue<-tissues_names$full_name[i]
path<-tissues_names$file_name[i]
gtex<-"/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"
gtex_tissue<- fread(paste0(gtex,path))



#----------------------
k=5
for(k in 1:length(unique(study_df$sample_id_rep))){
  print(k)
  sam<-unique(study_df$sample_id_rep)[k]
  
  ase_df<-study_df %>% filter(sample_id_rep==sam) 
  colnames(ase_df)[2]<-"pos"

xx<-sapply(strsplit(sam,"-"), function(xx){ paste0(xx[1],"-",xx[2]) })

gtex_tissue_1<-as_tibble(gtex_tissue[which(gtex_tissue$SUBJECT_ID %in% xx),])
colnames(gtex_tissue_1)[1:2]<-c("chr","pos")

join_ase<- gtex_tissue_1 %>% inner_join(ase_df)

join_ase<- join_ase %>% mutate(ratio_ref=log2(REF_COUNT) - log2(ref),
                                 mean_ref=(log2(REF_COUNT) + log2(ref))/2,
                                 ratio_alt=log2(ALT_COUNT) - log2(alt),
                                 mean_alt=(log2(ALT_COUNT) + log2(alt))/2,
                                 ratio_total=log2(TOTAL_COUNT) - log2(coverage),
                                 mean_total=(log2(TOTAL_COUNT) + log2(coverage))/2)


pdf(file=paste0("~/plot/ASE/GTEx_vs_recount/",sam,".pdf"), width = 10, height = 6)

p= ggplot(join_ase) + 
  geom_histogram(aes(REF_RATIO), fill="darkblue", alpha=0.3)+
  geom_histogram(aes(ref_ratio), fill="salmon", alpha=0.3)+
  labs(title="Same SNPs that are in Recount3: blue is GTEx, pink is Recount",
       subtitle=sam)+
  annotate("text", x= 0.7, y=3500, label= paste0("total gtex=" , nrow(gtex_tissue_1), ", total recount=", nrow(ase_df)))+
  annotate("text", x= 0.7, y=3000, label= paste0("gtex ref-ratio=" , round(median(join_ase$REF_RATIO,na.rm=T),3), ", recount ref-ratio=",round(median(join_ase$ref_ratio,na.rm=T),3)))
print(p)
  

p= ggplot(join_ase)+
  geom_point(aes(x=log10(ALT_COUNT), y=log10(REF_COUNT)), color="darkblue", alpha=0.3)+
  geom_point(aes(x=log10(alt), y=log10(ref)), color="salmon", alpha=0.3)
print(p)

p=ggplot(join_ase,aes(y=ratio_total,x=mean_total))+
  geom_point(alpha=0.4)+
  geom_smooth()+
  geom_hline(yintercept = 0,linetype = "dashed", alpha=0.5, color="red")
print(p)

p=ggplot(join_ase,aes(y=ratio_alt,x=mean_alt))+
  geom_point(alpha=0.4)+
  geom_smooth()+
  geom_hline(yintercept = 0,linetype = "dashed", alpha=0.5, color="red")
print(p)

p=ggplot(join_ase,aes(y=ratio_ref,x=mean_ref))+
  geom_point(alpha=0.4)+
  geom_smooth()+
  geom_hline(yintercept = 0,linetype = "dashed", alpha=0.5, color="red")
print(p)

p=ggplot(join_ase)+
  geom_point(aes(x=REF_RATIO,y=ref_ratio), alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="salmon")
print(p)

dev.off()
}


