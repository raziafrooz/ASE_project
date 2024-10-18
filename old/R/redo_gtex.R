setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(ggplot2)


geno_met<-read.csv("data/GTEx_geno_metadata.csv")
for(i in 10:length(geno_met$study)){
  study<-geno_met$study[i]
  print(study)
  
study_df<-as_tibble(readRDS(geno_met$allGenotypesOutput[geno_met$study==study]))%>% 
  filter(pred_genotype==2, coverage>=8)
study_df<-study_df %>% mutate(alt=round((sqrt(2^((2*S)-M))-1),0),ref=round((((2^M)*(alt+1))-1),0))
study_df$ref_ratio<- study_df$ref/study_df$coverage


saveRDS(study_df, file=paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rds") )
}
#-----------------------------------------------------






#------------------------------------------------------
#Get the true ASE hits from gtex 
#------------------------------------------------------
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



#-----------------------------------------------------
#Compare true GTEx to Reocunt3
#-----------------------------------------------------
tissues<-c("Lung","Stomach", "Liver", "Thyroid","Brain_Frontal_Cortex")
for(ss in 1:length(tissues)){
  study<-tissues[ss]
  print(study)
  study_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rds"))
  
  i<- which(tissues_names$full_name==study)
  print(tissues_names$full_name[i])
  abb<-tissues_names$abb[i]
  tissue<-tissues_names$full_name[i]
  path<-tissues_names$file_name[i]
  gtex<-"/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"
  gtex_tissue<- fread(paste0(gtex,path))
  
  
  
  l<-length(unique(study_df$sample_id_rep))
com_plot2<-data.frame(tissue=rep(study,l ),
                     sample=unique(study_df$sample_id_rep), 
                     recount3=rep(NA,l ),
                     gtex=rep(NA,l ))
for(k in 1:nrow(com_plot2)){
  sam<-com_plot2$sample[k]

  ase_df<-study_df %>% filter(sample_id_rep==sam) 
  colnames(ase_df)[2]<-"pos"

  xx<-sapply(strsplit(sam,"-"), function(xx){ paste0(xx[1],"-",xx[2]) })

  gtex_tissue_1<-as_tibble(gtex_tissue[which(gtex_tissue$SUBJECT_ID %in% xx),])
  colnames(gtex_tissue_1)[1:2]<-c("chr","pos")
  
  join_ase<- gtex_tissue_1 %>% inner_join(ase_df)
  
  com_plot2$recount3[k]<-median(join_ase$ref_ratio,na.rm=T)
  com_plot2$gtex[k]<-median(join_ase$REF_RATIO,na.rm=T)
  
  
}
if(ss==1){
  com_plot<-com_plot2
}else{
  com_plot<-rbind(com_plot,com_plot2)
}
}


  
  pdf(file="~/plot/ASE/recountVSgtex.pdf", width = 10, height = 6)
  ggplot(com_plot, aes(y=recount3, x=gtex, color=tissue))+
    geom_point(alpha=0.5)+
    ylim(c(0.47,0.52))+
    geom_vline(xintercept = 0.5,linetype = "dashed", alpha=0.5, color="red")+
    geom_hline(yintercept = 0.5,linetype = "dashed", alpha=0.5, color="red")+
    annotate("text", x = 0.505, y = 0.51, label = paste0(sum(com_plot$recount3==0.5, na.rm=T)," samples\nhave r-ratio=0.5 in Recount"))+
    labs(title="head to head comparison between Recount3 and GTEx",
         subtitle=paste0("total samples = ",nrow(com_plot)))
  dev.off()

  
  
  
  
  
  
  
#------------------------------------------------------------------
#plot
#------------------------------------------------------------------
  join_ase<- join_ase %>% mutate(ratio_ref=log2(REF_COUNT) - log2(ref),
                                   mean_ref=(log2(REF_COUNT) + log2(ref))/2,
                                   ratio_alt=log2(ALT_COUNT) - log2(alt),
                                   mean_alt=(log2(ALT_COUNT) + log2(alt))/2,
                                   ratio_total=log2(TOTAL_COUNT) - log2(coverage),
                                   mean_total=(log2(TOTAL_COUNT) + log2(coverage))/2)
  
  
  pdf(file="~/plot/ASE/GTEx_vs_recount/test.pdf", width = 10, height = 6)
  
  p= ggplot(join_ase) + 
    geom_histogram(aes(REF_RATIO), fill="darkblue", alpha=0.3)+
    geom_histogram(aes(ref_ratio), fill="salmon", alpha=0.3)
  print(p)
  p= ggplot(join_ase)+
    geom_point(aes(x=log10(ALT_COUNT), y=log10(REF_COUNT)), color="darkblue", alpha=0.3)+
    geom_point(aes(x=log10(alt), y=log10(ref)), color="salmon", alpha=0.3)
  print(p)
  
  p1 <- ggplot(join_ase) + 
    geom_point(aes(x=log10(ALT_COUNT), y=log10(REF_COUNT)), color="lightblue")
  p2 <- ggplot(join_ase) + 
    geom_point(aes(x=log10(alt), y=log10(ref)), color="salmon")
  
  p=plot_grid(p1, p2, labels = c('GTEx', 'Recount'))
  print(p)
  
  
  p1 <- ggplot(join_ase) + 
    geom_boxplot(aes(x=sample_id, y=REF_RATIO), color="lightblue")
  p2 <- ggplot(join_ase) + 
    geom_boxplot(aes(x=sample_id, y=ref_ratio), color="salmon")
  
  p=plot_grid(p1, p2, labels = c('GTEx', 'Recount'))
  print(p)
  
  p=ggplot(join_ase,aes(y=ratio_total,x=mean_total))+
    geom_point(alpha=0.4)+
    geom_smooth()
  print(p)
  
  p=ggplot(join_ase,aes(y=ratio_alt,x=mean_alt))+
    geom_point(alpha=0.4)+
    geom_smooth()
  print(p)
  
  p=ggplot(join_ase,aes(y=ratio_ref,x=mean_ref))+
    geom_point(alpha=0.4)+
    geom_smooth()
  print(p)
  
  p=ggplot(join_ase)+
    geom_point(aes(x=REF_RATIO,y=ref_ratio), alpha=0.4)+
    geom_abline(intercept = 0, slope = 1, color="salmon")
  print(p)
  
  dev.off()







