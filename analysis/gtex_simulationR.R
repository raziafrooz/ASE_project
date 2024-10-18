setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
#------------------------------------------------------
#Get the true ASE hits from gtex 
#------------------------------------------------------
#met<-read.csv("data/GTEx_metadata.csv")

tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("~/ASE-data/data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

tissues_names$file_name<-paste0("/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/", tissues_names$file_name)


#-----------------------------------------------------
#Compare true GTEx to Reocunt3
#-----------------------------------------------------
plotFile<-"/dcs07/hansen/data/recount_ASE/data/gtex_simulation.csv.gz"
if(!file.exists(plotFile)){
  tissues<-tissues_names$full_name
  all_mapp<-c()
  for(ss in 4:length(tissues)){
    study<-tissues[ss]
    print(study)
    
    gtex_tissue<- fread(tissues_names$file_name[tissues_names$full_name==study]) %>%
      select(CHR,POS,LOW_MAPABILITY,MAPPING_BIAS_SIM) %>% unique()
    colnames(gtex_tissue)[1:2]<- c("chr", "start")

    all_mapp<-rbind(all_mapp,gtex_tissue)
    
    all_mapp<-unique(all_mapp)
    
  }
  fwrite(all_mapp,plotFile)
  }
    
    