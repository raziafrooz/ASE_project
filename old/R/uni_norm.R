setwd("~/ASE/")
library(tidyverse)
library(ggplot2)
library(ggridges)
library(data.table)

single<-readRDS("data/single.rds")


#---------------------------------------------------
data<-single
data$uni_norm<-NA
for(k in 1:nrow(data)){
  sample_id<-data$external_id[k]
  study<-data$study[k]
  print(k)
  
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
    if(exists("ase_all")){
      ase_df<-ase_all
      rm(ase_all)
    } else{
      print("ase found")}
    ase_df<- ase_df %>% mutate(ratio=log10(alt) - log10(ref),
                               mean=(log10(alt) + log10(ref))/2)
    
    
    
    ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)
    
    q_line<-ase_df %>% group_by(mean_20tile) %>%
      reframe(median_ratio=median(ratio),
              enframe(quantile(ratio, c(0.05,0.95)), "quantile", "ratio_q"),
              min=min(mean),max=max(mean)) %>%
      mutate(q=case_when(quantile== '5%' ~ "low",
                         quantile== '95%' ~ "high" ))
    
    
    # positive_invert<- q_line[which(q_line$q=="high"),] %>%
    #   mutate(dist= ratio_q - median_ratio,
    #          invert= median_ratio-dist)
    # 
    # line_segment<-ase_df %>% group_by(mean_20tile) %>%
    #   summarise(min=min(mean),max=max(mean))
    
    test_line_sample<-left_join(test_line,q_line[,-6])
    test_line_sample<-test_line_sample %>% filter(q=="high",
                                                  CI=="high-ci",
                                                  max >=1,
                                                  max <=2.5)
    
    test_line_sample$geu_mean[which(test_line_sample$q=="high")]<- test_line_sample$geu_mean[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
    test_line_sample$CI_val[which(test_line_sample$q=="high")]<- test_line_sample$CI_val[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
    
    
    test_line_sample<-test_line_sample %>% mutate(diff= ratio_q-CI_val) %>%  filter(diff>0)
    
    data$uni_norm[k]<- mean(test_line_sample$diff,na.rm=T)
    
    rm(ase_df)
  }}
