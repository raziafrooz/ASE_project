library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library("org.Hs.eg.db")
library(VGAM)

source("~/ASE_project/src/remove_problematic_SNPs.R")
temp_dir<-"/dcs07/hansen/data/recount_ASE/tmp/"

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
#all_ase<-fread("~/test/geno_error.csv.gz")
for(t_id in 4:12){
  tissue_name<-tissues[t_id]
  print(tissue_name)
  
  if(!file.exists(paste0("~/test/geno_error_",tissue_name,".csv.gz"))){
   
gtex_prediction<-fread(gtex_test$genotypedSamples[gtex_test$tissue==tissue_name])
samples<-unique(gtex_prediction$sample_id_rep)
rm_id<-which(str_sub(samples, end= -3) %in% bad_samples)
if(length(rm_id)>0){samples<-samples[-rm_id]}

print(paste0("number of sample:",length(samples)))
all_ase<-c()
for(row_id in 1:length(samples)){
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
        
        ase_df$adj_ref<-ase_df$coverage-ase_df$adj_alt
        
        ase_df<-ase_df %>%dplyr::select(chr, start,adj_alt,adj_ref, coverage,true_genotype,pred_genotype,sample_id_rep)
        
        all_ase<-rbind(all_ase,ase_df)   
     
    }

fwrite(all_ase,paste0("~/test/geno_error_",tissue_name,".csv.gz"))
}
}

all_ase<-c()
for(t_id in 4:12){
  tissue_name<-tissues[t_id]
  print(tissue_name)
  
  ase<-fread(paste0("~/test/geno_error_",tissue_name,".csv.gz"))
  
  all_ase<-rbind(all_ase,ase)
}

#Get Prior 
yy<-all_ase %>% filter(true_genotype==1,coverage>20)
median_ratio=round(median(yy$adj_ref/yy$coverage),3) #0.857



calc_weight<-all_ase %>% group_by(sample_id_rep) %>%  summarize(true_1=sum(true_genotype==1)/n())
weight2<-median(calc_weight$true_1)
weight1<-1-weight2

get_error_p<-function(ref_count,coverage,weight1=0.98,weight2=0.02){
  nom<-(dbinom(ref_count, size = coverage, prob = 0.5)*dbinom((coverage-ref_count), size = coverage, prob = 0.5)*weight1 )
  err<-(dbinom(ref_count, size = coverage, prob = median_ratio)*dbinom((coverage-ref_count), size = coverage, prob =abs( 1-median_ratio))*weight2)
  p_val_error<- nom/(nom+err)
  return(p_val_error)
}

all_ase$genotyping_conf<-sapply(1:nrow(all_ase), function(zz) get_error_p(all_ase$adj_ref[zz],all_ase$coverage[zz]))
library(scattermore)
all_ase[1:4,]
pdf(file="~/plot/ASE/geno_error_2.pdf", width = 10, height = 6)
ggplot(data=all_ase ,aes(y=log2(adj_alt), x=log2(adj_ref)))+
  geom_scattermore(alpha=0.8,pointsize=2 )+
  geom_abline(slope=1, color="red")+
  geom_scattermore(data=all_ase %>% filter(true_genotype==1),aes(y=log2(adj_alt), x=log2(adj_ref)),aalpha=0.8,pointsize=2, color="purple")+
  geom_scattermore(data=all_ase %>% filter(genotyping_conf<0.00001),aes(y=log2(adj_alt), x=log2(adj_ref)),alpha=0.8,pointsize=2, color="blue")+
  labs(title="233 samples in 4 tissues. blue is homo alt, purple is homo ref")

  
ggplot(data=all_ase %>% filter(true_genotype==2),aes(adj_ref/coverage))+
  geom_histogram(alpha=0.4)+
  geom_histogram(data=all_ase %>% filter(true_genotype==1),aes(adj_ref/coverage),fill="purple",alpha=0.4)+
  geom_histogram(data=all_ase %>% filter(true_genotype==3),aes(adj_ref/coverage),fill="blue",alpha=0.4)+
  labs(title="233 samples in 4 tissues. gray is true het/purple is true homo ref/ blue is true homo alt")

dev.off()

#================================
#ROC curve
vv=0.1
try$error_perc<-aa
rr<-all_ase

roc_df<-data.frame(cut_off=c(0.0000000001,0.000000001,0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001))
for(i in 1:nrow(roc_df)){
  vv<-roc_df$cut_off[i]
  print(vv)
rr$rmv<-"no"
rr$rmv[rr$genotyping_conf<vv]<-"yes"


roc_df$true_pos[i]<-sum(rr$true_genotype==1 & rr$rmv=="yes") /sum(rr$true_genotype==1)
roc_df$false_pos[i]<-sum(rr$true_genotype==2 & rr$rmv=="yes") /sum(rr$true_genotype==2)

}

pdf(file="~/plot/ASE/geno_error_roc2.pdf", width = 10, height = 6)
ggplot(roc_df, aes(x=false_pos,y=true_pos,group=1,label=cut_off))+
  geom_line()+
  geom_point()+geom_text(hjust=0, vjust=0)+
  labs(title="ROC curve for genotyping error cut-off")
dev.off()

#-------------------------------------------------------
# one_sam<-fread(paste0("~/test/geno_error_",tissue_name,".csv.gz"))
# i=2
# one_sam<-one_sam %>% filter(sample_id_rep==unique(one_sam$sample_id_rep)[i])
# 
# one_sam$genotyping_conf<-sapply(1:nrow(one_sam), function(zz) get_error_p(one_sam$adj_ref[zz],one_sam$coverage[zz],pi1,pi2))
# one_sam[1,]
# 
# table(one_sam$true_genotype,round(one_sam$genotyping_conf2,1))
# 
# pdf(file="~/plot/ASE/test4.pdf", width = 10, height = 6)
# conf=0.000000001
# ggplot(one_sam %>% filter(genotyping_conf>conf), aes(x=log2(adj_ref),y=log2(adj_alt)))+
#   geom_point()+
#   geom_point(data=one_sam %>% filter(genotyping_conf<=conf), aes(x=log2(adj_ref),y=log2(adj_alt)),color="red", alpha=0.5)+
#   geom_point(data=one_sam %>% filter(one_sam$genotyping_conf<=conf,true_genotype==1), aes(x=log2(adj_ref),y=log2(adj_alt)),color="green3", alpha=0.6)+
#   geom_point(data=one_sam %>% filter(one_sam$genotyping_conf>conf,true_genotype==1), aes(x=log2(adj_ref),y=log2(adj_alt)),color="purple", alpha=0.6)+
#   labs(title=paste0(sum(one_sam$genotyping_conf<=conf&one_sam$true_genotype!=1,na.rm=T), " SNPs are het but removed (red)"),
#        subtitle=paste0(unique(one_sam$sample_id_rep),": (green) genotype confidence is <=",conf,". Purple is remaining geno error"))
# 
# dev.off()
# 
# sum(one_sam$genotyping_conf<=conf&one_sam$true_genotype!=1,na.rm=T)
# one_sam %>% filter(true_genotype==1,coverage>10) %>% summarize(m=median(adj_ref/coverage))
# 
# 
# yy[1,]
# median(yy$adj_ref/yy$coverage)
# median(yy$adj_alt/yy$coverage)
# table(one_sam$true_genotype,round(one_sam$genotyping_conf,1))
# 
