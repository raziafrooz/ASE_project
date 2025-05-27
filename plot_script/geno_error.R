library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library("org.Hs.eg.db")
library(VGAM)
library(cowplot)

theme_set(theme_cowplot())

source("~/ASE_project/src/remove_problematic_SNPs.R")

gtex_predict<-read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/metadata/all_gtex_metadata.csv")
gtex_test<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/GTEx_testing.csv")
gtex_uni_norm<-fread("/dcs07/hansen/data/recount_ASE/data/gtex_uniNorm.csv")
bad_samples<-gtex_uni_norm$sample_id[which(gtex_uni_norm$uni_norm_mean>0.08)]

#============================================================
#Read in wasp indv:
path_id<-"/dcl01/hansen/data/arazi/ASE/dbGap/phe000039.v1.GTEx_v8_ASE_WASP.expression-matrixfmt-ase.c1/GTEx_Analysis_v8_ASE_WASP_counts_by_subject/"
tissues_names<-as.data.frame(list.files(path =path_id))
colnames(tissues_names)<-"file_name"
tissues_names$indv<-gsub("\\..*","",tissues_names$file_name)
tissues_names$file_name<-paste0(path_id, tissues_names$file_name)
#============================================================
tissues<-unique(gtex_test$tissue)
t_id=37

  tissue_name<-tissues[t_id]
  print(tissue_name)

  
  gtex_prediction<-fread(gtex_test$genotypedSamples[gtex_test$tissue==tissue_name])
  samples<-unique(gtex_prediction$sample_id_rep)
  rm_id<-which(str_sub(samples, end= -3) %in% bad_samples)
  if(length(rm_id)>0){samples<-samples[-rm_id]}
  
  print(paste0("number of sample:",length(samples)))
  all_ase<-c()
 # for(row_id in 1:length(samples)){
  row_id=1
    print(row_id)
    sample_id<-samples[row_id]
    
    ase_df<-gtex_prediction%>% 
      filter(sample_id_rep== sample_id,pred_genotype==2, coverage>=8) %>% 
      mutate(ref_ratio=ref_count/coverage,
             alt_count=coverage-ref_count,
             ratio=log2(ref_count/alt_count),
             mean=(log2(ref_count)+log2(alt_count))/2)
    
    bigWig_path<-gtex_predict$total[which(gtex_predict$sample_id_rep==sample_id)][1]
    
    
    ase_df<-remove_problematic_SNPs(ase_df=ase_df,bigWig_path=bigWig_path)
    
    #Fix the counts by adjusting the MA plot to mean around 0 
    ratio_adj<-median(ase_df$ratio)
    adj<-ase_df$ratio-ratio_adj
    ase_df$adj_alt<-round((ase_df$coverage/((2^adj)+1)),0)
    
    stopifnot(sum(is.na(ase_df$adj_alt))==0)
    all_ase<-ase_df %>% mutate(adj_ref<-coverage-adj_alt)%>%
      dplyr::select(chr, start,adj_alt,adj_ref, coverage,true_genotype,pred_genotype,sample_id_rep)
 
    
    calc_weight<-all_ase %>% group_by(sample_id_rep) %>%  summarize(true_1=sum(true_genotype==1)/n())
    weight2<-median(calc_weight$true_1)
    weight1<-1-weight2
    
    get_error_p<-function(ref_count,coverage,median_ratio=0.86,weight1=0.98,weight2=0.02){
      nom<-(dbinom(ref_count, size = coverage, prob = 0.5)*dbinom((coverage-ref_count), size = coverage, prob = 0.5)*weight1 )
      err<-(dbinom(ref_count, size = coverage, prob = median_ratio)*dbinom((coverage-ref_count), size = coverage, prob =abs( 1-median_ratio))*weight2)
      p_val_error<- nom/(nom+err)
      return(p_val_error)}
    
    all_ase$genotyping_conf<-sapply(1:nrow(all_ase), function(zz) get_error_p(all_ase$adj_ref[zz],all_ase$coverage[zz]))
    library(scattermore)

p1<-ggplot(data=all_ase ,aes(y=log2(adj_alt), x=log2(adj_ref)))+
      geom_point(alpha=0.8)+
      geom_abline(slope=1, color="red")+
      geom_point(data=all_ase %>% filter(true_genotype==1),aes(y=log2(adj_alt), x=log2(adj_ref)),alpha=0.8, color="steelblue1")+
      geom_point(data=all_ase %>% filter(genotyping_conf<0.00001),aes(y=log2(adj_alt), x=log2(adj_ref)),alpha=0.8,color="magenta2")+
      labs(title=sample_id, y="Log2(Alt_count)",x="Log2(Ref_count)")+
  annotate("text", x=9.5, y=2.5, label= "True homozygous ref", color="magenta2") +
  annotate("text", x=9.5, y=2, label= "True homozygous alt", color="steelblue1")  
  
    
p2<-ggplot(data=all_ase %>% filter(true_genotype==2),aes(adj_ref/coverage))+
  geom_histogram(alpha=0.4)+
  geom_histogram(data=all_ase %>% filter(true_genotype==1),aes(adj_ref/coverage),fill="magenta2",alpha=0.4)+
  geom_histogram(data=all_ase %>% filter(true_genotype==3),aes(adj_ref/coverage),fill="steelblue1",alpha=0.4)+
  labs( y="Count",x="Reference ratio")


pdf(file="~/plot/ASE/geno_error.pdf", width = 10, height = 5)
plot_grid(p1,p2,labels = c('A', 'B'),align="hv")

dev.off() 


