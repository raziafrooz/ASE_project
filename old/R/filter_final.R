setwd("~/ASE/")
library(tidyverse)
library(ggplot2)
library(ggridges)
library(data.table)

single<-readRDS("data/single.rds")
paired<-readRDS("data/paired.rds")
cancer<-readRDS("data/cancer_annot.rds")
potential_single_cell<-readRDS("data/potential_single_cell.rds")





#--------------------------------------------------------
#Single end
#--------------------------------------------------------

scRNA<-which(single$external_id %in% potential_single_cell)

single<-single[-scRNA,]
low_count<-single %>% filter(read75>10)
nhet<- low_count %>% filter(nHet>500)
single_noncancer<-nhet[-which(nhet$sample_acc %in% cancer[,1]),]
fold_c<- single_noncancer %>% filter(ref_ratio>0.4, ref_ratio<0.6 )
good_single<- fold_c %>%  filter(uni_norm<0.185)




#--------------------------------------------------------
#Paired end
#--------------------------------------------------------
scRNA<-which(paired$external_id %in% potential_single_cell)

paired<-paired[-scRNA,]
low_count<-paired %>% filter(read75>10)
nhet<- low_count %>% filter(nHet>1000)
paired_noncancer<-nhet[-which(nhet$sample_acc %in% cancer[,1]),]
fold_c<- paired_noncancer %>% filter(ref_ratio>0.4, ref_ratio<0.6 )
good_paired<- fold_c %>%  filter(uni_norm<0.185)


good_ASE<-bind_rows(good_paired,good_single)
saveRDS(good_ASE, "data/good_ASE.rds")



pdf(file="~/plot/ASE/test/good_ASE.pdf", width = 10, height = 6)
ggplot(good_ASE, aes(ref_ratio, fill=library_layout))+
  geom_histogram(alpha=0.5)+
  labs(title= "Final samples that will be used for ASE",
       subtitle=paste0("#paired_end= ", nrow(good_paired), ", #single_end= ", nrow(good_single) ))
dev.off()



