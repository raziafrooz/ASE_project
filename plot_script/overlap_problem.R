library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
theme_set(theme_cowplot())

gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)

#===
samples=c("GTEX-1HBPH-0726-SM-ARZNA", "GTEX-N7MS-2526-SM-2D7W3", "GTEX-S4UY-0006-SM-3K2A7")


all_ase_df<-c()
for(i in samples){
row_id<-which(gtex_metadata$sample_id_rep == i)
sample_id<-gtex_metadata$sample_id_rep[row_id]

ase_df<-fread(gtex_metadata$genotypedSamples[which(gtex_metadata$sample_id_rep==sample_id)][1]) %>% 
  filter(pred_genotype==2, coverage>=8) %>% 
  mutate(ref_ratio=ref_count/coverage,
         alt_count=coverage-ref_count,
         ratio=log2(ref_count/alt_count),
         mean=(log2(ref_count)+log2(alt_count))/2)

ase_df$sample<-i
all_ase_df<-rbind(all_ase_df,ase_df)

}

all_ase_df$sample<-factor(all_ase_df$sample, levels=c("GTEX-S4UY-0006-SM-3K2A7","GTEX-1HBPH-0726-SM-ARZNA", "GTEX-N7MS-2526-SM-2D7W3"))


all_ase_df %>% group_by(sample) %>% summarize(median(ref_ratio))

pdf(file="~/plot/ASE/overlap_problem.pdf", width = 15, height = 6)
ggplot(all_ase_df,aes(log2(ref_count),log2(alt_count)) )+
  geom_point(alpha=0.4)+
  geom_abline(slope=1, linetype=2,color="red")+
  facet_grid(.~sample)+
  labs(y="log2(Alt)", x="log2(Ref)")
dev.off()





