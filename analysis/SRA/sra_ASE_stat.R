
library(tidyverse)
library(data.table)
# metadata<-read.csv("~/ASE/snakemake_ASE/ASE_metadata.csv")
# metadata[1,]
# metadata$ASE_path<-paste0("/dcs07/hansen/data/recount_ASE/output_SRA/",metadata$experiment_acc,"_ase.csv.gz")
# 
# fwrite(metadata,"/dcs07/hansen/data/recount_ASE/metadata/ASE_metadata.csv")
metadata<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/ASE_metadata.csv")

if(!file.exists("/dcs07/hansen/data/recount_ASE/data/sra_ASE_stat.csv.gz")){
stat_df<-c()
for (i in 1:nrow(metadata)){
  print(i)
  exp_id <- metadata$experiment_acc[i]
  
  ase_df<-fread(metadata$ASE_path[i])
  if(nrow(ase_df)>0){
  stat_df$experiment_acc[i]<-exp_id
  stat_df$nHet[i]<-nrow(ase_df)
  stat_df$ref_ratio[i]<-median(ase_df$ref_ratio)
  stat_df$q25[i]<-as.numeric(quantile(ase_df$coverage)[2])
  stat_df$q75[i]<-as.numeric(quantile(ase_df$coverage)[4])
  }else{
    stat_df$experiment_acc[i]<-exp_id
    stat_df$nHet[i]<-nrow(ase_df)
    stat_df$ref_ratio[i]<-0
    stat_df$q25[i]<-0
    stat_df$q75[i]<-0
    
  }
}


df<-as.data.frame(stat_df)
fwrite(df,"/dcs07/hansen/data/recount_ASE/data/sra_ASE_stat.csv.gz")
}else{
  df<-fread("/dcs07/hansen/data/recount_ASE/data/sra_ASE_stat.csv.gz")
}


