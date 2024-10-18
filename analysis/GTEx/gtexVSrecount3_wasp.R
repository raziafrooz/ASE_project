
setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(MetBrewer)

tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("~/ASE-data/data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

tissues_names$file_name<-paste0("/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/", tissues_names$file_name)


#-----------------------------------------------------
#Compare true GTEx to Reocunt3
#-----------------------------------------------------
#Get gtex genotype metadata:
gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)
plotFile<-"/dcs07/hansen/data/recount_ASE/data/gtexVSrecount_wasp.csv.gz"
if(!file.exists(plotFile)){
  tissues<-tissues_names$full_name
  
  for(ss in 1:length(tissues)){
    study<-tissues[ss]
    print(study)
    
    gtex_tissue<- fread(tissues_names$file_name[tissues_names$full_name==study][1])
    colnames(gtex_tissue)[1:2]<- c("chr", "start")
    
    l<-length(unique(gtex_tissue$SAMPLE_ID))
    com_plot2<-data.frame(tissue=rep(study,l ),
                          sample=unique(gtex_tissue$SAMPLE_ID), 
                          recount3=rep(NA,l ),
                          gtex=rep(NA,l ))
    for(k in 1:nrow(com_plot2)){
      sam<-com_plot2$sample[k]
      
      
      ase_df<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==sam][1]) %>%
        filter(pred_genotype==2, coverage>=8) %>% 
        mutate(ref_ratio=ref_count/coverage)
      
      gtex_tissue_sam<-gtex_tissue %>% filter(SAMPLE_ID==sam) 
      
      
      
      join_ase<- inner_join(ase_df,gtex_tissue_sam)
      
      com_plot2$recount3[k]<-median(join_ase$ref_ratio,na.rm=T)
      com_plot2$gtex[k]<-median(join_ase$REF_RATIO,na.rm=T)
      
      
    }
    if(ss==1){
      com_plot<-com_plot2
    }else{
      com_plot<-rbind(com_plot,com_plot2)
    }
  }
  fwrite(com_plot,plotFile)
}else{
  com_plot<-fread(plotFile)
}



pdf(file="~/plot/ASE/recountVSgtex_wasp.pdf", width = 10, height = 6)

pp=ggplot(com_plot, aes(y=recount3, x=gtex, color=tissue))+
  geom_point(alpha=0.3)+
  ylim(c(0.47,0.52))+
  geom_vline(xintercept = 0.5,linetype = "dashed", alpha=0.5, color="red")+
  geom_hline(yintercept = 0.5,linetype = "dashed", alpha=0.5, color="red")+
  annotate("text", x = 0.505, y = 0.51, label = paste0(sum(round(com_plot$recount3,2)==0.5, na.rm=T)," samples\nhave r-ratio=0.5 in Recount"))+
  labs(title="head to head comparison between Recount3 and GTEx (WASP QC passed)",
       subtitle=paste0("total samples = ",sum(!is.na(com_plot$recount3))))+
  theme(legend.position="none")

print(pp)


ggplot(com_plot, aes(recount3)) + 
  geom_density()+
  labs(title="Recount3 ref_ratio density plot (WASP QC passed)")


dev.off()








