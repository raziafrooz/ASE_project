setwd("~/ASE/")
library(tidyverse)
library(ggplot2)
library(ggridges)
library(data.table)
library(rmarkdown)
library(tinytex)

render("R/filter_steps.Rmd")
render("R/filter_steps_single.Rmd")

paired<-readRDS("data/paired.rds")
battle<-read_tsv("data/battle.tsv")
test_line<-readRDS("data/geuvadis_quantile.rds")
cancer<-readRDS("data/cancer_annot.rds")
potential_single_cell<-readRDS("data/potential_single_cell.rds")



paired$frag_mode<-paired$bc_frag.mode_length
paired$overlap<- paired$bc_frag.mode_length-(paired$avg_len*2)



paired<-paired[-which(paired$study %in% potential_single_cell),]



seq_mean=seq(0,4.6,by=0.1)

source("scr/MA_plot.R")

#------------------------------------------------------------------------------------
#filter based on ref_ratio and nHet
#------------------------------------------------------------------------------------
# paired$nHet<-NA
# paired$ref_ratio<-NA
#  paired$read75<-NA
# 
# for(k in 1:nrow(paired)){
#   sample_id<-paired$external_id[k]
#   study<-paired$study[k]
#   print(k)
#   if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
#   load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
#   if(exists("ase_all")){
#     ase_df<-ase_all
#     rm(ase_all)
#   } else{
#     print("ase found")}
# 
#     paired$nHet[k]<-nrow(ase_df)
#     paired$ref_ratio[k]<-median(ase_df$ref_ratio)
#     paired$read75[k]<-as.numeric(quantile(ase_df$total,0.75))
#     
#     rm(ase_df)
#   }}
#saveRDS(paired, "data/paired.rds")
#------------------------------------------------------------------------------------
paired_noncancer<-paired[-which(paired$sample_acc %in% cancer[,1]),]
low_count<-paired_noncancer %>% filter(read75>=10)
nhet<- low_count %>% 
  filter(nHet>1000)
#8780
nrow(paired)-nrow(paired_noncancer)
nrow(paired_noncancer)- nrow(low_count)
  nrow(low_count)-nrow(nhet)

data<-nhet
data_samp<-data[sample(nrow(data), 100),]
pdf(file="~/plot/ASE/test/samp_nhet.pdf", width = 10, height = 6)
for(k in 1:nrow(data_samp)){
  sample_id<-data_samp$external_id[k]
  study<-data_samp$study[k]
  uni_norm<-round(data_samp$uni_norm[k],3)
  #overlap<-data_samp$overlap[k]
  n_het<-data_samp$nHet[k]
  print(k)
  
  xx<-make_MA(sample_id,study,uni_norm,test_line)
  
  p1<-plot_MA(xx$ase_df, sample_id,study,uni_norm,xx$test_line_sample,xx$q_line)
  p1<-p1+#annotate("text", x = 4, y = 0.8, label = paste0("fragment overlap=",overlap) )+
    annotate("text", x = 4, y = 0.7, label = paste0("nHet=",n_het) )
  print(p1)
  
}
dev.off()
#------------------------------------------------------------------------------------

filter<-nhet %>% 
  filter(ref_ratio>0.4, ref_ratio<0.6 )
nrow(nhet)-nrow(filter)
-nrow(good)

data_samp<-data[sample(nrow(data), 100),]
pdf(file="~/plot/ASE/test/samp_ratio.pdf", width = 10, height = 6)
for(k in 1:nrow(data_samp)){
  sample_id<-data_samp$external_id[k]
  study<-data_samp$study[k]
  uni_norm<-round(data_samp$uni_norm[k],3)
  #overlap<-data_samp$overlap[k]
  n_het<-data_samp$nHet[k]
  print(k)
  
  xx<-make_MA(sample_id,study,uni_norm,test_line)
  
  p1<-plot_MA(xx$ase_df, sample_id,study,uni_norm,xx$test_line_sample,xx$q_line)
  p1<-p1+#annotate("text", x = 4, y = 0.8, label = paste0("fragment overlap=",overlap) )+
    annotate("text", x = 4, y = 0.7, label = paste0("nHet=",n_het) )
  print(p1)
  
}
dev.off()

#------------------------------------------------------------------------------------
good<-filter %>% filter(uni_norm<0.185)
dim(filter)-dim(good)

data_samp<-good[sample(nrow(good), 100),]
pdf(file="~/plot/ASE/test/samp_good.pdf", width = 10, height = 6)
for(k in 1:nrow(data_samp)){
  sample_id<-data_samp$external_id[k]
  study<-data_samp$study[k]
  uni_norm<-round(data_samp$uni_norm[k],3)
  overlap<-data_samp$overlap[k]
  #ref_ratio<-round(data_samp$ref_ratio[k],2)
  n_het<-data_samp$nHet[k]
  print(k)
  
  xx<-make_MA(sample_id,study,uni_norm,test_line)
  
  p1<-plot_MA(xx$ase_df, sample_id,study,uni_norm,xx$test_line_sample,xx$q_line)
  p1<-p1+annotate("text", x = 4, y = 0.8, label = paste0("fragment overlap=",overlap) )+
    annotate("text", x = 4, y = 0.7, label = paste0("nHet=",n_het) )
  print(p1)
  
}
dev.off()
#------------------------------------------------------------------------------------


