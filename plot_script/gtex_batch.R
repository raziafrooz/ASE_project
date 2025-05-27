library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)

theme_set(theme_cowplot())
theme1= theme(axis.text=element_text(size = 10), axis.title=element_text(size = 20))+
  theme_half_open() +
  background_grid()
  
  
qc_df<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/gtex_qc_metadata.csv.gz")
qc_df<-qc_df %>% mutate(overlap= (star.average_input_read_length)-bc_frag.mode_length)
qc_df$SAMPLE_ID<-str_sub(qc_df$external_id, end= -3)

gtex_recount<-read.csv("/dcs07/hansen/data/recount_ASE/data/gtex_recountPipeline.csv.gz")
qc_df$ref_ratio<-gtex_recount$ref_ratio[match(qc_df$external_id,gtex_recount$sample_id)]


pdf(file="~/plot/ASE/gtex_batch.pdf", width = 9, height = 5)
ggplot(qc_df,aes(x=overlap,y=ref_ratio,color=SMGEBTCHT))+
  geom_point(alpha=0.5)+theme1+
  theme(legend.position=c(.05,.15))+
  labs(title="GTEx samples processed through Recount3 pipeline",y="Reference Ratio", x="Overlap")
dev.off()

