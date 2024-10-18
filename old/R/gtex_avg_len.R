
library(recount3)
setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(ggplot2)

study<-"Lung"

#--------------------------------------------
#get true gtex:
met<-read.csv("data/GTEx_metadata.csv")
tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]


#--------------------------------------------
studies<- met$tissue[1:20]
i=1

for(k in 1:length(studies)){
  print(k)
  study<-studies[k]
study_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rds"))
study_df$sample_id<-sapply(strsplit(study_df$sample_id_rep,"-"), function(xx){ paste0(xx[1],"-",xx[2]) })


i<- which(tissues_names$full_name==study)
print(tissues_names$full_name[i])
abb<-tissues_names$abb[i]
tissue<-tissues_names$full_name[i]
path<-tissues_names$file_name[i]
gtex<-"/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"
gtex_tissue<- fread(paste0(gtex,path))




colnames(gtex_tissue)[c(1,2,7)]<-c("chr","start","sample_id")
ase_joined<-inner_join(gtex_tissue,study_df,by=c("chr","start","sample_id"))
ase_joined_sum<-ase_joined %>% group_by(sample_id_rep) %>% summarize(gtex=median(REF_RATIO), recount=median(ref_ratio))
ase_joined_sum$tissue<-study

if(study==studies[1]){
  plot<-ase_joined_sum
} else{
  plot<-rbind(plot,ase_joined_sum)
}
}


human_projects <- available_projects()
qc_meta<-subset(human_projects, file_source == "gtex" & project_type == "data_sources")
plot$frag_len<-NA

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


plot$frag_len[id2]<-xx$bc_frag.mean_length[id]

xx <-utils::read.delim(file_retrieve(url[4], verbose = FALSE))

id<-match(plot$sample_id_rep,xx$external_id)
id<-id[!is.na(id)]
id2<-match(xx$external_id,plot$sample_id_rep)
id2<-id2[!is.na(id2)]

plot$avg_len[id2]<-xx$avg_len[id]
plot$num_bases[id2]<-xx$X.bases[id]

}
#saveRDS(plot, "data/fragment_plot.rds")
plot$seq_overlap<- plot$frag_len-(plot$avg_len*2)
# plot[which(plot$recount<0.4),]
# xx[1:4,]
# xx[xx$external_id%in%c("GTEX-S32W-0326-SM-5S2T4.1","GTEX-VUSG-0926-SM-5S2SV.1","GTEX-T6MN-0826-SM-5S2SK.1"),]
quantile(plot$seq_overlap)
plot_final[which(plot_final$ref_ratio<0.38),]
plot[1,]
plot_final<-plot %>% pivot_longer(c(recount,gtex), names_to = "dataset", values_to = "ref_ratio")

pdf(file="~/plot/ASE/gtex_seq_len_more.pdf", width = 10, height = 6)

ggplot(plot_final)+
  geom_violin(aes(x=factor(avg_len),y=ref_ratio, fill=dataset))+
  geom_jitter(data= plot_final[plot_final$dataset=="recount",],
              aes(x=factor(avg_len),y=ref_ratio, color=tissue), alpha=0.5,width = 0.15)+
  labs(title="average sequencing length for 2,100 gtex samples",
       subtitle="for the same SNPs, ref ratio was calculated in gtex and reocunt3 pipeline")
dev.off()



plot_final<-plot_final %>% mutate(seq_overlap_cut= cut_number(seq_overlap,10))

plot_final<-plot %>% mutate(ref_ratio_cut=cut_number(recount,2)) %>% 
  pivot_longer(c(recount,gtex), names_to = "dataset", values_to = "ref_ratio") %>% 
  mutate(fg_cut=cut_number(frag_len,10),
         seq_overlap_cut=cut_number(seq_overlap,10),
         num_bases_cut=cut_number(num_bases,10))
  



library(ggridges)
plot_final<-plot_final %>%   mutate(fg_cut=cut_number(frag_len,10),
                                    seq_overlap_cut=cut_number(seq_overlap,10))
plot2<-plot %>% mutate(seq_overlap_cut=cut_number(seq_overlap,10))

pdf(file="~/plot/ASE/seq_overlap.pdf", width = 10, height = 6)

ggplot(plot_final)+
  geom_violin(aes(x=seq_overlap_cut,y=ref_ratio, fill=dataset))+
  geom_jitter(data= plot_final[plot_final$dataset=="recount",],
              aes(x=seq_overlap_cut,y=ref_ratio, color=tissue), alpha=0.5,width = 0.15)+
  labs(title="frag_len-(avg_len*2) for 1,877 gtex samples (20tissues)",
       subtitle="for the same SNPs, ref ratio of reocunt3 pipeline is plotted as points")

ggplot(plot_final %>% filter(ref_ratio>=0.45))+
  geom_violin(aes(x=seq_overlap_cut,y=ref_ratio, fill=dataset))+
  geom_jitter(data= plot_final[plot_final$dataset=="recount",] %>% filter(ref_ratio>=0.45),
              aes(x=seq_overlap_cut,y=ref_ratio, color=tissue), alpha=0.5,width = 0.15)+
  labs(title="NO OUTLIER:frag_len-(avg_len*2) for 1,877 gtex samples (20tissues)",
       subtitle="for the same SNPs, ref ratio of reocunt3 pipeline is plotted as points")


ggplot(plot_final[plot_final$dataset=="recount",])+
  geom_density_ridges(aes(x=ref_ratio,y=seq_overlap_cut, fill=dataset),alpha=0.6)+
  xlim(c(0.46,0.52))+
  labs(title="NO OUTLIER:frag_len-(avg_len*2) for 1,877 gtex samples (20tissues)")

ggplot(plot2)+
  geom_point(aes(gtex,recount, color=seq_overlap_cut), alpha=0.7)+
  facet_wrap(vars(seq_overlap_cut))

dev.off()

plot_final$group<-"high"
plot_final$group[plot_final$ref_ratio=>0.4 & plot_final$ref_ratio <= 0.46 ]<-"med1"
plot_final$group[plot_final$ref_ratio>0.46 & plot_final$ref_ratio <= 0.47 ]<-"med2"
plot_final$group[plot_final$ref_ratio>0.47 & plot_final$ref_ratio <= 0.495 ]<-"med3"
plot_final$group[plot_final$ref_ratio<0.4]<-"low"

pdf(file="~/plot/ASE/test.pdf", width = 10, height = 6)

ggplot(plot_final[plot_final$dataset=="recount",])+
  geom_density_ridges(aes(x=num_bases/1e+6,y=group, fill=dataset),alpha=0.5)

dev.off()

as.data.frame(plot[plot$recount>0.46& plot$recount<0.475,])[1:3,]
as.data.frame(plot[plot$recount<0.46,])[1:3,]
as.data.frame(plot[plot$recount==0.5,])[1:3,]

as.data.frame(plot[which(plot$mapped_multi_loci==max(plot$mapped_multi_loci)),])

plot$too_many_mismatches<-NA
plot$unmap_other<-NA
plot$mapped_multi_loci<-NA
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
  
  #%_of_reads_mapped_to_multiple_loci: Number of reads mapped to multiple loci divided by number of input reads
  col_num<-which(str_detect(colnames(xx),"reads_mapped_to_multiple_loci"))[1]
  plot$mapped_multi_loci[id2]<-xx[,col_num][id]
  
  #%_of_reads_unmapped:_too_many_mismatches: Number of reads where best alignment has more mismatches than max allowed number of mismatches divided by number of input reads
  col_num<-which(str_detect(colnames(xx),"_of_reads_unmapped._too_many_mismatches"))[1]
  plot$too_many_mismatches[id2]<-xx[,col_num][id]
  
  
  #%_of_reads_unmapped:_other: Reads are unmapped due to no acceptable seed/windows divided by number of input reads
  col_num<-which(str_detect(colnames(xx),"_of_reads_unmapped._other"))[1]
  plot$unmap_other[id2]<-xx[,col_num][id]
  
}



plot$monorail_all_base<-NA
plot$monorail_unique<-NA
plot$morotail_percent<-NA
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
  
  #%_of_reads_mapped_to_multiple_loci: Number of reads mapped to multiple loci divided by number of input reads
  col_num<-which(str_detect(colnames(xx),"gene_fc_count_unique.total"))[1]
  plot$monorail_all_base[id2]<-xx[,col_num][id]
  
  #%_of_reads_unmapped:_too_many_mismatches: Number of reads where best alignment has more mismatches than max allowed number of mismatches divided by number of input reads
  col_num<-which(str_detect(colnames(xx),"gene_fc.unique_"))[1]
  plot$monorail_unique[id2]<-xx[,col_num][id]
  
  
  #%_of_reads_unmapped:_other: Reads are unmapped due to no acceptable seed/windows divided by number of input reads
  col_num<-which(str_detect(colnames(xx),"gene_fc_count_unique.assigned"))[1]
  plot$morotail_percent[id2]<-xx[,col_num][id]
  
}

as.data.frame(plot[1,])

plot2<-plot %>% mutate(monorail_unique_cut=cut_number(monorail_unique,10),
                       monorail_all_base_cut=cut_number(monorail_all_base,10),
                       morotail_percent_cut=cut_number(morotail_percent,10),
                       recount_cut=cut_number(recount,2))



pdf(file="~/plot/ASE/test.pdf", width = 10, height = 6)

ggplot(plot2 %>% filter(recount>=0.45))+
  geom_violin(aes(x=monorail_unique_cut,y=recount))+
  geom_jitter(aes(x=monorail_unique_cut,y=recount), alpha=0.5,width = 0.15)
  

ggplot(plot2 %>% filter(recount>=0.45))+
  geom_violin(aes(x=monorail_all_base_cut,y=recount))+
  geom_jitter(aes(x=monorail_all_base_cut,y=recount), alpha=0.5,width = 0.15)

ggplot(plot2 %>% filter(recount>=0.45))+
  geom_violin(aes(x=morotail_percent_cut,y=recount))+
  geom_jitter(aes(x=morotail_percent_cut,y=recount), alpha=0.5,width = 0.15)

dev.off()



