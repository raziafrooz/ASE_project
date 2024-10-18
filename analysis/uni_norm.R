#geuvadis ERP001942 
library(quantreg)
library(data.table)
library(tidyverse)
setwd("~/ASE/")
sra_met<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_SRA.csv")
geu<-sra_met[which(sra_met$study=="ERP001942"),]
bad_samples<-c("ERR188390","ERR204881","ERR204916","ERR204940",
               "ERR204843","ERR204897","ERR205016","ERR188328")

ase_metadata<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/ASE_metadata.csv")
ase_metadata<-ase_metadata[which(ase_metadata$sample_id %in% geu$sample_id),]
if(!file.exist("geuvadis_quantile_new.rds")){
#Train based on Geuvadis:
seq_mean=seq(0,4.6,by=0.1)
q_line_geu<-c()
for(k in 1:nrow(ase_metadata)){
  sample_id<-ase_metadata$sample_id[k]
  study<-ase_metadata$study[k]
  print(k)
  if(sample_id %in% bad_samples)
    { print("bad sample")}else{
    
    ase_df<-fread(ase_metadata$ASE_path[k]) %>% 
      mutate(ratio=log10(alt_count) - log10(ref_count),
             mean=(log10(alt_count) + log10(ref_count))/2)
    
    
    
    ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)
    
    
    
    q_line<-ase_df %>% group_by(mean_20tile) %>%
      reframe(median_ratio=median(ratio),
                enframe(quantile(ratio, c(0.05,0.95)), "quantile", "ratio_q"),
                min=min(mean),max=max(mean)) %>% 
      mutate(q=case_when(quantile== '5%' ~ "low",
                         quantile== '95%' ~ "high" ))
    q_line$sample<-sample_id
    

      q_line_geu<-rbind(q_line_geu,q_line)
    }

}


test_line<-q_line_geu %>%
  group_by(mean_20tile,q) %>%
  reframe(geu_mean=mean(ratio_q),
            enframe(quantile(ratio_q, c(0.05,0.95)), "CI", "CI_val")) %>%
  mutate(CI=case_when(CI== '5%' ~ "low_ci",
                     CI== '95%' ~ "high-ci" ))
test_line<-test_line[!is.na(test_line$mean_20tile),]
test_line$max<-sort(c(seq_mean[-c(1,2,3,4,5)],seq_mean[-c(1,2,3,4,5)],seq_mean[-c(1,2,3,4,5)],seq_mean[-c(1,2,3,4,5)]))
saveRDS(test_line, "geuvadis_quantile_new.rds")

}else{
  test_line<-readRDS("geuvadis_quantile_new.rds")
}
