library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)

theme_set(theme_cowplot())
theme1= theme(axis.text=element_text(size = 10), axis.title=element_text(size = 20))


qc_df<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/gtex_qc_metadata.csv.gz")
qc_df<-qc_df %>% mutate(overlap= (star.average_input_read_length)-bc_frag.mode_length)
qc_df$SAMPLE_ID<-str_sub(qc_df$external_id, end= -3)

gtex_recount<-read.csv("/dcs07/hansen/data/recount_ASE/data/gtex_recountPipeline.csv.gz")
qc_df$ref_ratio<-gtex_recount$ref_ratio[match(qc_df$external_id,gtex_recount$sample_id)]

gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)

id<-which(qc_df$ref_ratio<0.4)[1]
sample_bad<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==qc_df$SAMPLE_ID[id]])
sample_bad<-sample_bad %>% filter(coverage>=8,pred_genotype==2) %>% mutate(alt_count=coverage-ref_count)
sample_bad$sample<-qc_df$SAMPLE_ID[id]


id<-which(qc_df$ref_ratio==0.5)[2]
sample_good<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==qc_df$SAMPLE_ID[id]])
sample_good<-sample_good %>% filter(coverage>=8,pred_genotype==2) %>% mutate(alt_count=coverage-ref_count)
sample_good$sample<-qc_df$SAMPLE_ID[id]

sam_plot<-rbind(sample_good,sample_bad)

pdf(file="~/plot/ASE/contamination.pdf", width = 8.5, height = 5)
ggplot(sam_plot,aes(x=log2(ref_count),y=log2(alt_count)))+
  geom_point(alpha=0.5)+
  geom_abline(slope=1,color="red")+
  facet_wrap(vars(sample))+theme1
dev.off()
#--------------------------------------------------------------------
#Overlap
theme1= theme(axis.text=element_text(size = 10), axis.title=element_text(size = 20))+
  theme_half_open() +
  background_grid()

id<-which(qc_df$overlap< -80 & qc_df$overlap> -100)[50]
sample_bad<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==qc_df$SAMPLE_ID[id]])
sample_bad<-sample_bad %>% filter(coverage>=8,pred_genotype==2) %>%
  mutate(alt_count=coverage-ref_count,
         ratio=log2(ref_count/alt_count),
         mean=(log2(ref_count)+log2(alt_count))/2)
sample_bad$sample<-qc_df$SAMPLE_ID[id]


id<-which(qc_df$overlap> 50 & qc_df$overlap<90)[1]
sample_good<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==qc_df$SAMPLE_ID[id]])
sample_good<-sample_good %>% filter(coverage>=8,pred_genotype==2) %>% 
  mutate(alt_count=coverage-ref_count,
         ratio=log2(ref_count/alt_count),
         mean=(log2(ref_count)+log2(alt_count))/2)
sample_good$sample<-qc_df$SAMPLE_ID[id]

sam_plot<-rbind(sample_good,sample_bad)

id<-which(qc_df$overlap== 0 )[1]
sample_noOV<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==qc_df$SAMPLE_ID[id]])
sample_noOV<-sample_noOV %>% filter(coverage>=8,pred_genotype==2) %>%
  mutate(alt_count=coverage-ref_count,
         ratio=log2(ref_count/alt_count),
         mean=(log2(ref_count)+log2(alt_count))/2)
sample_noOV$sample<-qc_df$SAMPLE_ID[id]
sam_plot<-rbind(sam_plot,sample_noOV)

sam_plot$sample<-factor(sam_plot$sample, levels=c("GTEX-S4UY-0006-SM-3K2A7",
                                                  "GTEX-1HBPH-0726-SM-ARZNA",
                                                  "GTEX-N7MS-2526-SM-2D7W3"))

pdf(file="~/plot/ASE/overlap.pdf", width = 10, height = 5)
ggplot(sam_plot,aes(x=mean,y=ratio))+
  geom_point(alpha=0.3)+
  geom_hline(yintercept=0,color="red")+
  facet_wrap(sample~.)+
  theme1+
  geom_smooth()+
  labs(x="A",y="M")
dev.off()








