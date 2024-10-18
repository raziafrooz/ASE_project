setwd("~/ASE/")
library(tidyverse)
library(ggplot2)
library(scattermore)

#########################
#Plot data analysis: choose 3 tissues 
#########################
tissues<-c("Liver","Stomach","Lung")
for(tissue in tissues){
  print(tissue)
ase_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", tissue, ".rds"))
true_gtex<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/true_gtex_", tissue, ".rds"))
colnames(true_gtex)[1:2]<-c("chr","pos")
#Inner join to only select the locations that are same between the two dataset
ase_join<-inner_join(true_gtex,ase_df, by=c("chr", "pos", "sample_id"))
ase_join$tissue<-tissue
if(tissue == tissues[1]){
  ase_plot<-ase_join
} else {
  ase_plot<-rbind(ase_plot, ase_join)
}
}

ase_plot$ref_ratio<- ase_plot$ref/ase_plot$total

plot <- ase_plot %>% 
  group_by(tissue, sample_id) %>%
  summarize(Recount_median=median(ref_ratio), GTEx_median= median(REF_RATIO))

pdf(file="~/plot/ASE/gtexVSrecount_3tissue.pdf", width = 10, height = 4)
ggplot(plot %>% filter(Recount_median > 0.35), aes(y= Recount_median, x=GTEx_median, color=tissue))+
  geom_point(alpha=0.5)+
  labs(title="Comparison of ref_ratio between GTEx and Recount3 in 3 tissues")
  #geom_hline(yintercept = 0.5)+
  #geom_vline(xintercept = 0.5)
dev.off()

plot<-ase_plot %>% rowwise() %>% mutate(mine_allele=min(alt,ref),MIN_ALLELE= min(REF_COUNT,ALT_COUNT))

pdf(file="~/plot/ASE/gtexVSrecount_mono_allelic.pdf", width = 10, height = 4)

  p=ggplot(plot)+
    geom_scattermore(aes(x=total, y=mine_allele), alpha=0.4)+
    xlim(0,2000)+ylim(0,1000)
  print(p)
  
  p=ggplot(plot)+
    geom_scattermore(aes(x=TOTAL_COUNT, y=MIN_ALLELE), alpha=0.4)+
    xlim(0,2000)+ylim(0,1000)
  print(p)
  
dev.off()




