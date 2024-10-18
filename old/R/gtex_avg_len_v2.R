
library(recount3)
setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(ggplot2)
library(ggridges)

gtex<-readRDS("data/fragment_plot.rds")
gtex



human_projects <- available_projects()
qc_meta<-subset(human_projects, file_source == "gtex" & project_type == "data_sources")
gtex$frag_len_mod<-NA

for(ss in qc_meta$project){
  print(ss)
  url<-locate_url(
    ss,
    "data_sources/gtex",
    type = "metadata")
  
  xx <-utils::read.delim(file_retrieve(url[3], verbose = FALSE))
  
  id<-match(gtex$sample_id_rep,xx$external_id)
  id<-id[!is.na(id)]
  id2<-match(xx$external_id,gtex$sample_id_rep)
  id2<-id2[!is.na(id2)]
  
  
  gtex$frag_len_mod[id2]<-xx$bc_frag.mode_length[id]
  
}
  
gtex$overlap<-gtex$frag_len_mod-(2*gtex$avg_len)
gtex<-gtex %>% mutate(overlap_cut=cut_number(overlap,4),
                frag_mod_cut= cut_number(frag_len_mod,4),
                ref_ratio_cut=cut_number(recount,2))

pdf(file="~/plot/ASE/gtex_overlap.pdf", width = 10, height = 6)

p= ggplot(gtex)+
  geom_violin(aes(x=overlap_cut,y=recount))+
  geom_jitter(aes(x=overlap_cut,y=recount), alpha=0.5,width = 0.15)+
  labs(title="1,877 GTEx samples from 20 tissues: overlap is [frag_len-(2*seq_len)] ")
print(p)

p= ggplot(gtex)+
  geom_point(aes(x=overlap,y=recount))+
  labs(title="1,877 GTEx samples from 20 tissues: overlap is [frag_len-(2*seq_len)] ")
print(p)


p= ggplot(gtex)+
  geom_violin(aes(x=frag_mod_cut,y=recount))+
  geom_jitter(aes(x=frag_mod_cut,y=recount), alpha=0.5,width = 0.15)+
  labs(title="1,877 GTEx samples from 20 tissues: mod fragment length")
print(p)
dev.off()








  