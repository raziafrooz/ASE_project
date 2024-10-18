setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(ggplot2)
library(MetBrewer)

cc_fill <-scale_fill_manual(values=met.brewer("Hokusai1"))
cc_color <-scale_color_manual(values=met.brewer("Hokusai1"))

#Get gtex genotype metadata:
gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
#met_testing<-read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/GTEx_testing.csv")
#met_training<-read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/GTEx_training.csv")


plotFile<-"/dcs07/hansen/data/recount_ASE/data/gtex_recountPipeline.csv.gz"
if(!file.exists(plotFile)){
plot<-c()
for(i in 1:nrow(gtex_metadata)){
  
  tissue<-gtex_metadata$study[i]
  sample<-gtex_metadata$sample_id[i]
  print(tissue)
  
  study_df<-as_tibble(fread(gtex_metadata$genotypedSamples[i]))%>% 
    filter(pred_genotype==2, coverage>=8) 
  
  plot2<-data.frame(sample_id=sample,
             tissue,
             ref_ratio=median(study_df$ref_count/study_df$coverage))
  
    plot<-rbind(plot,plot2)
  }
  
fwrite(plot,plotFile)
}else{
  plot<-fread(plotFile)
}

plot<-plot %>% mutate(indv_id= paste0(str_split_i(sample_id,"-",1),"-",str_split_i(sample_id,"-",2)))
bad_indv<-unique(plot$indv_id[which(plot$ref_ratio<0.45)])


pdf(file="~/plot/ASE/GTEx_recount_ase.pdf", width = 10, height = 4)
ggplot(plot,aes(y=ref_ratio, x=tissue))+
  geom_point(alpha=0.3)+ 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))+
  geom_hline(yintercept = 0.5,linetype = "dashed", alpha=0.5, color="red")+
  labs(title=paste0("samples in all gtex tissues: total = ", nrow(plot)),
       subtitle=paste0( sum(round(plot$ref_ratio,2)==0.5)," samples have exactly ref-ratio = 0.5")) +
  cc_fill+
  cc_color

ggplot(plot[plot$ref_ratio>0.45,],aes(y=ref_ratio, x=tissue))+
  geom_point(alpha=0.3)+ 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))+
  geom_hline(yintercept = 0.5,linetype = "dashed", alpha=0.5, color="red")+
  labs(title=paste0("samples in all gtex tissues: total = ", nrow(plot)),
       subtitle=paste0( sum(round(plot$ref_ratio,2)==0.5)," samples have exactly ref-ratio = 0.5")) +
  cc_fill+
  cc_color

ggplot(plot[which(plot$indv_id %in% bad_indv),],aes(y=ref_ratio, x=tissue, color=indv_id))+
  geom_jitter()+ 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))+
  geom_hline(yintercept = 0.5,linetype = "dashed", alpha=0.5, color="red")+
  labs(title="Outliers are bad replicates of few samples in different tissues")
  # cc_fill+
  # cc_color
dev.off()




