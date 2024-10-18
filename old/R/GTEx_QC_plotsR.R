
library(recount3)
setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(ggplot2)


geno_met<-read.csv("data/GTEx_geno_metadata.csv")
human_projects <- available_projects()
qc_meta<-subset(human_projects, file_source == "gtex" & project_type == "data_sources")
#plot<-readRDS("data/all_gtex_plot.rds")
plot<-readRDS("data/fragment_plot.rds")

plot$frag_len<-NA
ss<-qc_meta$project[1]
for(ss in qc_meta$project){
  print(ss)
  url<-locate_url(
  ss,
  "data_sources/gtex",
  type = "metadata")

xx <-utils::read.delim(file_retrieve(url[3], verbose = FALSE))

id<-match(plot$sample_id_rep,xx$external_id)
id<-id[!is.na(id)]
id2<-match(xx$external_id,plot$sample_id_rep)
id2<-id2[!is.na(id2)]



lung<-plot[id2,]
xx$ref_ratio<-NA
xx$ref_ratio[id]<-lung$recount
xx<-xx[!is.na(xx$ref_ratio),]


if(ss == qc_meta$project[1]){
  all_qc<-xx}else{
all_qc<-rbind(all_qc,xx)
  }
}
#saveRDS(all_qc, file= "data/gtex_qc_plot.rds" )
all_qc<-readRDS(all_qc)

pdf(file="~/plot/ASE/test.pdf", width = 10, height = 6)
for(name_col in colnames(all_qc)[4:113]){
  if( sum(duplicated(quantile(all_qc[,name_col]))) == 0   ){
    name_col <- colnames(all_qc)[15]
  print(name_col)
  plot2<- all_qc %>% mutate(cut=cut_number(.data[[name_col]],20))


  
 p= ggplot(plot2)+
    geom_violin(aes(x=cut,y=ref_ratio))+
    geom_jitter(aes(x=cut,y=ref_ratio), alpha=0.5,width = 0.15)+
    labs(title=name_col,
         x=paste0(name_col, "_cut"))
  print(p)
  }else{ 
    print("cannot make the plot")
  }
}
dev.off()


colnames(all_qc)[]
dim(all_qc)
sum(all_qc$bc_frag.kallisto_mean_length>=195 & all_qc$bc_frag.kallisto_mean_length<=250 )

#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
plot$avg_len<-NA
plot$frag_mod<-NA

ss<-qc_meta$project[1]
for(ss in qc_meta$project){
  print(ss)
  url<-locate_url(
    ss,
    "data_sources/gtex",
    type = "metadata")
  
  xx <-utils::read.delim(file_retrieve(url[4], verbose = FALSE))
  
  id<-match(plot$sample_id_rep,xx$external_id)
  id<-id[!is.na(id)]
  id2<-match(xx$external_id,plot$sample_id_rep)
  id2<-id2[!is.na(id2)]
  
  
  plot$avg_len[id2]<-xx$avg_len[id]
  
  
  xx <-utils::read.delim(file_retrieve(url[3], verbose = FALSE))
  
  id<-match(plot$sample_id_rep,xx$external_id)
  id<-id[!is.na(id)]
  id2<-match(xx$external_id,plot$sample_id_rep)
  id2<-id2[!is.na(id2)]
  
  
  plot$frag_mod[id2]<-xx$bc_frag.mode_length[id]
}

plot$mode_overlap<- plot$frag_mod-(plot$avg_len*2)



pdf(file="~/plot/ASE/gtex_mode_frag.pdf", width = 10, height = 6)


plot3<-plot %>%  mutate( ref_ratio_cut=cut_number(ref_ratio,2))
p= ggplot(plot3)+
  geom_density_ridges(aes(x=frag_mod,y=ref_ratio_cut), alpha=0.5)+
  labs(title="gtex using mode.frag")
print(p)


plot3<-plot %>%  mutate( mode_overlap_cut=cut_number(mode_overlap,7),
                         ref_ratio_cut=cut_number(ref_ratio,2))





p= ggplot(plot3)+
  geom_violin(aes(x=mode_overlap_cut,y=ref_ratio))+
  geom_jitter(aes(x=mode_overlap_cut,y=ref_ratio), alpha=0.5,width = 0.15)+
  labs(title="Using mode frag to calculate (mode.frag - avg_len*2)")
print(p)



dev.off()



