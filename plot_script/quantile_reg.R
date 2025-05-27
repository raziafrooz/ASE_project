library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
theme_set(theme_cowplot())
source("~/ASE_project/src/uniNorm_function.R")

test_line<-readRDS("~/ASE-data/geuvadis_quantile_new.rds")
gtex_uni<-fread(("/dcs07/hansen/data/recount_ASE/data/gtex_uniNorm.csv"))
gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)


plot_list<-list()
for(sample in c("GTEX-1117F-0226-SM-5GZZ7","K-562-SM-43V9A")){
k=which(gtex_uni$sample_id==sample)



indiv<-gtex_metadata %>% filter(sample_id_rep==gtex_uni$sample_id[k])

ase_df<-fread(indiv$genotypedSamples) %>%
  filter(pred_genotype==2, coverage>=8) %>% 
  mutate(ref_ratio=ref_count/coverage,
         alt_count=coverage-ref_count)



xx<-make_MA(ase_df,test_line)
study<-"GTEx"
sample_id=gtex_uni$sample_id[k]
uni_norm=gtex_uni$uni_norm_mean[k]

plot_list[[sample]]<-plot_MA(xx$ase_df, sample_id,study,uni_norm,xx$test_line_sample,xx$q_line)

}


pdf(file="~/plot/ASE/quantile_reg.pdf", width = 10.5, height = 5)

plot_grid(plotlist=plot_list,labels = c('A', 'B'),align="hv")


dev.off()






