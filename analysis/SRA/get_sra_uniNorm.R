library(data.table)
library(tidyverse)
source("~/ASE/src/uniNorm_function.R")

metadata<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/ASE_metadata.csv")
test_line<-readRDS("~/ASE/geuvadis_quantile_new.rds")
final_tissue<-fread("/dcs07/hansen/data/recount_ASE/data/toUse_tissue_primaryCell.csv")
#final_tissue<-fread("/dcs07/hansen/data/recount_ASE/data/toUse_no_annot.csv")
#paired<-final_tissue %>%   sample_n(30)
uni_norm_df<-c()
#pdf(file="~/plot/ASE/test.pdf", width = 10, height = 6)

for(i in 1:nrow(final_tissue)){
  print(i)
  exp_id<-final_tissue$experiment_acc[i]
  sample_id<-final_tissue$sample_accession[i]
xx<-make_MA(exp_id, test_line)

uni_norm<- xx$test_line_sample %>%
  filter(q=="high",CI=="high-ci",max >=1,max <=2.5) %>% 
  mutate(diff= ratio_q-CI_val) %>% 
  summarize(uni_norm_mean=mean(diff,na.rm=T))


uni_norm_df<-rbind(uni_norm_df,uni_norm)

#pp<-plot_MA(xx$ase_df, exp_id, sample_id,uni_norm$mm, xx$test_line_sample, xx$q_line)

#print(pp)
}
#dev.off()
uni_norm_df$experiment_acc<-final_tissue$experiment_acc

if(!file.exists("/dcs07/hansen/data/recount_ASE/data/ase_uniNorm.csv")){
fwrite(xx,"/dcs07/hansen/data/recount_ASE/data/ase_uniNorm.csv")
}
# 
# if(!file.exists("/dcs07/hansen/data/recount_ASE/data/ase_uniNorm_noAnno.csv")){
#   fwrite(xx,"/dcs07/hansen/data/recount_ASE/data/ase_uniNorm_noAnno.csv")
# }

