
library(recount3)
setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(ggplot2)

#----------------------
#ASE from GTEx pipeline:

tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

tissues_names$file_name<-paste0("/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/", tissues_names$file_name)

#--------------------

qc_df<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/gtex_qc_metadata.csv.gz")
gtex_recount<-read.csv("/dcs07/hansen/data/recount_ASE/data/gtex_recountPipeline.csv.gz")
qc_df$ref_ratio<-gtex_recount$ref_ratio[match(qc_df$external_id,gtex_recount$sample_id)]


plot_df<-qc_df


pdf(file="~/plot/ASE/gtex_qc_plots.pdf", width = 10, height = 6)
for(name_col in colnames(qc_df)[73:198]){
  
    print(name_col)

    pp=ggplot(qc_df,aes(x=.data[[name_col]], y= ref_ratio, color = study))+
      geom_point()+ theme(legend.position="none")+
      labs(title=name_col)

    print(pp)
  }
dev.off()

qc_df<-qc_df %>% mutate(overlap= (star.average_input_read_length)-bc_frag.mode_length)

pdf(file="~/plot/ASE/gtex_qc_plots_overlap.pdf", width = 10, height = 6)
  
  pp=ggplot(qc_df,aes(x=overlap, y= ref_ratio, color = study))+
    geom_point()+ theme(legend.position="none")+
    labs(title="Ovelap is (star.average_input_read_length-bc_frag.mode_length)")
  
  print(pp)
  
  
  pp=ggplot(qc_df,aes(y=overlap, x= SMGEBTCHT, color = study))+
    geom_jitter()+ theme(legend.position="none")+
    labs(title="Ovelap is (star.average_input_read_length-bc_frag.mode_length)")
  
  print(pp)
dev.off()

#------------------------------------------------------------------------------------------------
# compare some gtex samples to recount individualy 
#------------------------------------------------------------------------------------------------
com_plot<-read.csv("/dcs07/hansen/data/recount_ASE/data/gtexVSrecount.csv.gz")
gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)
com_plot$recount3<-round(com_plot$recount3,3)
com_plot$gtex<-round(com_plot$gtex,3)


make_plot<-function(recount_higher){
study<-unique(recount_higher$tissue)
print(study)

gtex_tissue<- fread(tissues_names$file_name[tissues_names$full_name==study])

colnames(gtex_tissue)[1:2]<- c("chr", "start")
sam<-recount_higher$sample
  
#Get recount genotype
  ase_df<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==sam][1]) %>%
    filter(pred_genotype==2, coverage>=8) %>% 
    mutate(ref_ratio=ref_count/coverage,
           alt_count=coverage-ref_count)
  #Get Gtex
  gtex_tissue_sam<-gtex_tissue %>% filter(SAMPLE_ID==sam) 
  
  join_ase<- inner_join(ase_df,gtex_tissue_sam) %>% 
    select(coverage,ref_count,alt_count,
           ref_ratio,REF_COUNT,ALT_COUNT,
           TOTAL_COUNT,REF_RATIO) %>% 
    mutate(ratio_total=log2(TOTAL_COUNT) - log2(coverage),
          mean_total=(log2(TOTAL_COUNT) + log2(coverage))/2,
          ratio_ref=log2(REF_COUNT) - log2(ref_count),
          mean_ref=(log2(REF_COUNT) + log2(ref_count))/2,
          ratio_alt=log2(ALT_COUNT) - log2(alt_count),
          mean_alt=(log2(ALT_COUNT) + log2(alt_count))/2)
  

  return(join_ase)
  }


recount_higher<-com_plot[which(com_plot$recount3<com_plot$gtex)[1],]
join_ase<-make_plot(recount_higher)

  pdf(file="~/plot/ASE/recountISbad.pdf", width = 10, height = 6)
  sam<-recount_higher$sample
  pp=ggplot(join_ase,aes(x=ref_ratio, y= REF_RATIO))+
    geom_point()+ theme(legend.position="none")+
    labs(title="REF_RATIO")+
    geom_smooth()+
    geom_abline(intercept = 0, slope = 1, color="red")+
    labs(title=paste0(sam,", nSNP= ", nrow(join_ase)),
         subtitle=paste0("gtex ref_ratio = " ,round(median(join_ase$REF_RATIO),3),
                         ", recount ref_ratio = " ,round(median(join_ase$ref_ratio),3), ", overlap=",qc_ss$overlap))
  
  print(pp)
  
  pp=ggplot(join_ase,aes(x=mean_total, y= ratio_total))+
    geom_point()+ theme(legend.position="none")+
    labs(title="Total")+
    geom_smooth()+
    geom_hline(yintercept = 0, color="red", linetype=2)+
    labs(title=paste0(sam,", nSNP= ", nrow(join_ase)),
         subtitle=paste0("gtex ref_ratio = " ,round(median(join_ase$REF_RATIO),3),
                         ", recount ref_ratio = " ,round(median(join_ase$ref_ratio),3), ", overlap=",qc_ss$overlap))
  
  
  print(pp)
  
  pp=ggplot(join_ase,aes(x=mean_ref, y= ratio_ref))+
    geom_point()+ theme(legend.position="none")+
    labs(title="Ref")+
    geom_smooth()+
    geom_hline(yintercept = 0, color="red", linetype=2)+
    labs(title=paste0(sam,", nSNP= ", nrow(join_ase)),
         subtitle=paste0("gtex ref_ratio = " ,round(median(join_ase$REF_RATIO),3), 
                         ", recount ref_ratio = " ,round(median(join_ase$ref_ratio),3), ", overlap=",qc_ss$overlap))
  
  
  print(pp)
  
  pp=ggplot(join_ase,aes(x=mean_alt, y= ratio_alt))+
    geom_point()+ theme(legend.position="none")+
    labs(title="Alt")+
    geom_smooth()+
    geom_hline(yintercept = 0, color="red", linetype=2)+
    labs(title=paste0(sam,", nSNP= ", nrow(join_ase)),
         subtitle=paste0("gtex ref_ratio = " ,round(median(join_ase$REF_RATIO),3), 
                         ", recount ref_ratio = " ,round(median(join_ase$ref_ratio),3), ", overlap=",qc_ss$overlap))
  
  
  print(pp)
  dev.off()

  pdf(file="~/plot/ASE/gtex_indv_ase2.pdf", width = 10, height = 6)
  
  pp=ggplot(join_ase,aes(x=log2(REF_COUNT), y= log2(ALT_COUNT)))+
    geom_point()+ theme(legend.position="none")+
    labs(title="Alt")+
    geom_smooth()+
    geom_abline(intercept = 0, slope = 1, color="red")+
    labs(title=paste0("Gtex counts for ",sam,", nSNP= ", nrow(join_ase)),
         subtitle=paste0("gtex ref_ratio = " ,round(median(join_ase$REF_RATIO),3), ", recount ref_ratio = " ,round(median(join_ase$ref_ratio),3)))
  
  print(pp)
  
  pp=ggplot(join_ase,aes(x=log2(ref_count), y= log2(alt_count)))+
    geom_point()+ theme(legend.position="none")+
    labs(title="Alt")+
    geom_smooth()+
    xlim(c(0,12))+
    geom_abline(intercept = 0, slope = 1, color="red")+
    labs(title=paste0("Recount counts for ",sam,", nSNP= ", nrow(join_ase)),
         subtitle=paste0("gtex ref_ratio = " ,round(median(join_ase$REF_RATIO),3), ", recount ref_ratio = " ,round(median(join_ase$ref_ratio),3)))
  
  print(pp)
  
  dev.off()
  
  
  
  
  
  #------------------------------

  recount_higher<-com_plot[which(com_plot$sample %in% qc_df$SAMPID[which(qc_df$overlap>50)]),][1,]
  #recount_higher<-com_plot[which(com_plot$recount3<com_plot$gtex)[1],]
  join_ase<-make_plot(recount_higher)
  
  pdf(file="~/plot/ASE/test_overlapPOSITIVE.pdf", width = 10, height = 6)
  
  sam<-recount_higher$sample
  qc_ss<-qc_df[which(qc_df$SAMPID==sam),]
  
  pp=ggplot(join_ase,aes(x=ref_ratio, y= REF_RATIO))+
    geom_point()+ theme(legend.position="none")+
    labs(title="REF_RATIO")+
    geom_smooth()+
    geom_abline(intercept = 0, slope = 1, color="red")+
    labs(title=paste0(sam,", nSNP= ", nrow(join_ase)),
         subtitle=paste0("gtex ref_ratio = " ,round(median(join_ase$REF_RATIO),3),
                         ", recount ref_ratio = " ,round(median(join_ase$ref_ratio),3), ", overlap=",qc_ss$overlap ))
  
  print(pp)
  
  pp=ggplot(join_ase,aes(x=mean_total, y= ratio_total))+
    geom_point()+ theme(legend.position="none")+
    labs(title="Total")+
    geom_smooth()+
    geom_hline(yintercept = 0, color="red", linetype=2)+
    labs(title=paste0(sam,", nSNP= ", nrow(join_ase)),
         subtitle=paste0("gtex ref_ratio = " ,round(median(join_ase$REF_RATIO),3), 
                         ", recount ref_ratio = " ,round(median(join_ase$ref_ratio),3), ", overlap=",qc_ss$overlap))
  
  
  print(pp)
  
  pp=ggplot(join_ase,aes(x=mean_ref, y= ratio_ref))+
    geom_point()+ theme(legend.position="none")+
    labs(title="Ref")+
    geom_smooth()+
    geom_hline(yintercept = 0, color="red", linetype=2)+
    labs(title=paste0(sam,", nSNP= ", nrow(join_ase)),
         subtitle=paste0("gtex ref_ratio = " ,round(median(join_ase$REF_RATIO),3), 
                         ", recount ref_ratio = " ,round(median(join_ase$ref_ratio),3), ", overlap=",qc_ss$overlap))
  
  
  print(pp)
  
  pp=ggplot(join_ase,aes(x=mean_alt, y= ratio_alt))+
    geom_point()+ theme(legend.position="none")+
    labs(title="Alt")+
    geom_smooth()+
    geom_hline(yintercept = 0, color="red", linetype=2)+
    labs(title=paste0(sam,", nSNP= ", nrow(join_ase)),
         subtitle=paste0("gtex ref_ratio = " ,round(median(join_ase$REF_RATIO),3),
                         ", recount ref_ratio = " ,round(median(join_ase$ref_ratio),3), ", overlap=",qc_ss$overlap))
  
  
  print(pp)
  dev.off()
  