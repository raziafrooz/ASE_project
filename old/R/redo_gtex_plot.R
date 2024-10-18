setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(ggplot2)

geno_met<-read.csv("data/GTEx_geno_metadata.csv")
for(i in 1:length(geno_met$study)){
  study<-geno_met$study[i]
  print(study)
  study_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rds") )
  
  plot2<-study_df %>% group_by(sample_id_rep) %>% summarize(ref_ratio=median(ref_ratio))
  plot2$tissue<-study
  
  if(i==1){
    plot<-plot2
  }else{
    plot<-rbind(plot,plot2)
  }
  
}

saveRDS(plot, file = "data/all_gtex_plot.rds")

pdf(file="~/plot/ASE/allGtex.pdf", width = 10, height = 4)
ggplot(plot,aes(y=ref_ratio, x=tissue))+
  geom_point(alpha=0.3)+ 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))+
  geom_hline(yintercept = 0.5,linetype = "dashed", alpha=0.5, color="red")+
  labs(title="test set samples in all gtex tissues (3,901 samples)",
       subtitle="2709 samples (69%) have ref-ratio = 0.5")

ggplot(plot,aes(y=ref_ratio, x=tissue))+
  geom_point(alpha=0.3)+ 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))+
  geom_hline(yintercept = 0.5,linetype = "dashed", alpha=0.5, color="red")+
  ylim(c(0.45,0.55))+
  labs(title="test set samples in all gtex tissues (3,901 samples)",
       subtitle="2709 samples (69%) have ref-ratio = 0.5")
dev.off()

