library(data.table)
library(tidyverse)
source("~/ASE/src/uniNorm_function.R")

test_line<-readRDS("~/ASE/geuvadis_quantile_new.rds")

gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)


uni_norm_df<-c()
#pdf(file="~/plot/ASE/test.pdf", width = 10, height = 6)

for(i in 1:nrow(gtex_metadata)){
  print(i)
  exp_id<-unique(gtex_metadata$sample_id_rep)[i]

  ase_df<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==exp_id][1]) %>%
    filter(pred_genotype==2, coverage>=8) %>% 
    mutate(ref_ratio=ref_count/coverage,
           alt_count=coverage-ref_count)

  xx<-make_MA(ase_df, test_line)
  
  uni_norm<- xx$test_line_sample %>%
    filter(q=="high",CI=="high-ci",max >=1,max <=2.5) %>% 
    mutate(diff= ratio_q-CI_val) %>% 
    summarize(uni_norm_mean=mean(diff,na.rm=T))
  
  
  uni_norm_df<-rbind(uni_norm_df,uni_norm)
  uni_norm_df$sample_id<-exp_id
  #pp<-plot_MA(xx$ase_df, exp_id, exp_id,uni_norm$uni_norm_mean, xx$test_line_sample, xx$q_line)
  
  #print(pp)
}
#dev.off()

if(!file.exists("/dcs07/hansen/data/recount_ASE/data/gtex_uniNorm.csv")){
     fwrite(uni_norm_df,"/dcs07/hansen/data/recount_ASE/data/gtex_uniNorm.csv")
   }
  
