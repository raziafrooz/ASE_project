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


id<-which(qc_df$overlap> 50 & qc_df$overlap<90)[2]
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

sam_plot$sample<-factor(sam_plot$sample, levels=c("GTEX-N7MS-2526-SM-2D7W3",
                                                  "GTEX-1HBPH-0726-SM-ARZNA",
                                                  "GTEX-S4UY-0006-SM-3K2A7"))

pdf(file="~/plot/ASE/overlap2.pdf", width = 10, height = 5)
ggplot(sam_plot,aes(x=mean,y=ratio))+
  geom_point(alpha=0.3)+
  geom_hline(yintercept=0,color="red")+
  facet_wrap(sample~.)+
  theme1+
  geom_smooth()+
  labs(x="A",y="M")
dev.off()

#---------------------------------------------------------------
#Fix overlap


source("~/ASE_project/src/remove_problematic_SNPs.R")

id<-which(qc_df$overlap> 50 & qc_df$overlap<90)[1]
sample_good<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==qc_df$SAMPLE_ID[id]])
sample_good<-sample_good %>% filter(coverage>=8,pred_genotype==2) %>% 
  mutate(alt_count=coverage-ref_count,
         ratio=log2(ref_count/alt_count),
         mean=(log2(ref_count)+log2(alt_count))/2)
sample_good$sample<-qc_df$SAMPLE_ID[id]

sample_good<-remove_problematic_SNPs(ase_df=sample_good)

#Fix the counts by adjusting the MA plot to mean around 0
ratio_adj<-median(sample_good$ratio)

  
adj<-sample_good$ratio-ratio_adj
sample_good$adj_alt<-round(sample_good$coverage/((2^adj)+1),0)
sample_good$adj_ref<-sample_good$coverage-sample_good$adj_alt

sample_good<-sample_good %>% mutate(ratio_corrected=log2(adj_ref/adj_alt),
                                    mean_corrected=(log2(adj_ref)+log2(adj_alt))/2)


before<-ggplot(sample_good,aes(x=mean,y=ratio))+
  geom_point(alpha=0.3)+geom_hline(yintercept=0,color="red")+theme1+geom_smooth()+labs(x="A",y="M")
after<-ggplot(sample_good,aes(x=mean_corrected,y=ratio_corrected))+
  geom_point(alpha=0.3)+ geom_hline(yintercept=0,color="red")+theme1+geom_smooth()+labs(x="A",y="M")

plot_row <- plot_grid(before,after,labels = "AUTO")

title <- ggdraw() + 
  draw_label(paste0(sample_good$sample," overlap correction") ,fontface = 'bold',x = 0,hjust = 0) +
  theme( plot.margin = margin(0, 0, 0, 7) )

pdf(file="~/plot/ASE/overlap_fix.pdf", width = 12, height = 5)
plot_grid(title,plot_row,ncol = 1, rel_heights = c(0.1, 1) )
dev.off()



