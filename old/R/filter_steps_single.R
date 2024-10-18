setwd("~/ASE/")
library(tidyverse)
library(ggplot2)
library(ggridges)
library(data.table)
source("scr/MA_plot.R")

single<-readRDS("data/single.rds")
battle<-read_tsv("data/battle.tsv")
test_line<-readRDS("data/geuvadis_quantile.rds")
cancer<-readRDS("data/cancer_annot.rds")
potential_single_cell<-readRDS("data/potential_single_cell.rds")


id<-which(single$sample_acc %in% cancer[,1])
length(id)
single<-single[-id,]
id2<-which(single$external_id %in% potential_single_cell)
length(id2)
single<-single[-id2,]

seq_mean=seq(0,4.6,by=0.1)
#------------------------------------------------------------------------------------
#filter based on ref_ratio and nHet
#------------------------------------------------------------------------------------
# single$nHet<-NA
# single$ref_ratio<-NA
#  single$read75<-NA
# 
# for(k in 1:nrow(single)){
#   sample_id<-single$external_id[k]
#   study<-single$study[k]
#   print(k)
#   if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
#   load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
#   if(exists("ase_all")){
#     ase_df<-ase_all
#     rm(ase_all)
#   } else{
#     print("ase found")}
# 
#     single$nHet[k]<-nrow(ase_df)
#     single$ref_ratio[k]<-median(ase_df$ref_ratio)
#     single$read75[k]<-as.numeric(quantile(ase_df$total,0.75))
# 
#     rm(ase_df)
#   }}
# saveRDS(single, "data/single.rds")
#------------------------------------------------------------------------------------
quantile(single$nHet,0.2,na.rm=T)

#low_count<-single_noncancer %>% filter(read75>=10)
nhet<- single %>% 
  filter(nHet>500)
#8780
nrow(single)-nrow(nhet)
nrow(single)- nrow(low_count)
nrow(low_count)-nrow(nhet)


data<-nhet
data_samp<-data[sample(nrow(data), 100),]
pdf(file="~/plot/ASE/test/single_samp_nhet.pdf", width = 10, height = 6)
for(k in 1:nrow(data_samp)){
  sample_id<-data_samp$external_id[k]
  study<-data_samp$study[k]
  #uni_norm<-round(data_samp$uni_norm[k],3)
  #overlap<-data_samp$overlap[k]
  n_het<-data_samp$nHet[k] 
  
  xx<-make_MA(sample_id,study,uni_norm,test_line)
  
  p1<-plot_MA(xx$ase_df, sample_id,study,uni_norm=1,xx$test_line_sample,xx$q_line)
  p1<-p1+#annotate("text", x = 4, y = 0.8, label = paste0("fragment overlap=",overlap) )+
    annotate("text", x = 4, y = 0.7, label = paste0("nHet=",n_het) )
  print(p1)
  
}
dev.off()

#------------------------------------------------------------------------------------
quantile(single$read75, na.rm=T)
low_count<-nhet %>% filter(read75>10)
filter<-low_count %>% 
  filter(ref_ratio>0.4, ref_ratio<0.6 )
nrow(nhet)-nrow(filter)
-nrow(good)

#------------------------------------------------------------------------------------

quantile(single$read75,na.rm=T)
#low_count<-single_noncancer %>% filter(read75>=10)
low_count<-nhet %>% filter(read75>10)
uni<-low_count %>%  filter(uni_norm>0.185)

data<-uni
data_samp<-data[sample(nrow(data), 100),]
pdf(file="~/plot/ASE/test/single_samp_uni.pdf", width = 10, height = 6)
for(k in 1:nrow(data_samp)){
  sample_id<-data_samp$external_id[k]
  study<-data_samp$study[k]
  #uni_norm<-round(data_samp$uni_norm[k],3)
  #overlap<-data_samp$overlap[k]
  n_het<-data_samp$nHet[k] 
  
  xx<-make_MA(sample_id,study,uni_norm,test_line)
  
  p1<-plot_MA(xx$ase_df, sample_id,study,uni_norm=1,xx$test_line_sample,xx$q_line)
  p1<-p1+#annotate("text", x = 4, y = 0.8, label = paste0("fragment overlap=",overlap) )+
    annotate("text", x = 4, y = 0.7, label = paste0("nHet=",n_het) )
  print(p1)
  
}
dev.off()


