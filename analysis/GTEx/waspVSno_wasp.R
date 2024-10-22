setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(MetBrewer)
gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)


tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("~/ASE-data/data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

tissues_names$file_name<-paste0("/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/", tissues_names$file_name)
wasp<-tissues_names


tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("~/ASE-data/data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

tissues_names$file_name<-paste0("/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/", tissues_names$file_name)
no_wasp<-tissues_names
rm(tissues_names)

plotFile<-"~/test/waspVSno_wasp.csv.gz"
com_plot<-c()
tissues<-wasp$full_name

if(!file.exists(plotFile)){
for(ss in 2:length(tissues)){
  study<-tissues[ss]
  print(study)
  
  wasp_1<- fread(wasp$file_name[wasp$full_name==study][1])
  colnames(wasp_1)[1:2]<- c("chr", "start")
  
  
  no_wasp_1<- fread(no_wasp$file_name[no_wasp$full_name==study][1])
  colnames(no_wasp_1)[1:2]<- c("chr", "start")
  
  
  
  
  
  join_ase<- inner_join(wasp_1,no_wasp_1, by=c("chr", "start","SAMPLE_ID"))
  rm(no_wasp_1)
  rm(wasp_1)
  
  com_plot2<- join_ase %>%
    group_by(SAMPLE_ID) %>% 
    summarize(wasp=median(REF_RATIO.x,na.rm=T),
              no_wasp=median(REF_RATIO.y,na.rm=T))
  
  
  com_plot<-rbind(com_plot,com_plot2)
  rm(com_plot2)
}
fwrite(com_plot,plotFile)
}else{
  
  com_plot<-fread(plotFile)
}



#===================================================================

quantile(plotFile_wasp$wap_recount)
quantile(plotFile_no_wasp$no_wap_recount)

plotFile_no_wasp<-fread("/dcs07/hansen/data/recount_ASE/data/gtexVSrecount.csv.gz")
plotFile_wasp<-fread("/dcs07/hansen/data/recount_ASE/data/gtexVSrecount_wasp.csv.gz")
colnames(plotFile_no_wasp)[3:4]<-c("no_wap_recount","no_wasp_gtex")
colnames(plotFile_wasp)[3:4]<-c("wap_recount","wasp_gtex")

plot_df<-inner_join(plotFile_wasp,plotFile_no_wasp)
plot_df<-plot_df %>%
  pivot_longer(!c(tissue,sample), names_to = "pipeline", values_to = "ref_ratio")


plot_df[which(plot_df$no_wap_recount>0.48 & plot_df$wap_recount<=0.48 &tissue==study),][1:5,]

pdf(file="~/plot/ASE/waspVSnowasp2.pdf", width = 10, height = 6)


ggplot(plot_df, aes(x=pipeline, y=round(ref_ratio,2)))+
  geom_boxplot()+
  labs(title="Boxplot of median ref ratio of samples in GTEx with and without WASP correction",
       subtitle="recount was subseted to common SNPs with GTEx in each condition")
dev.off()


#===================================================================
# 
# 
# ss=1
# tissues<-wasp$full_name
# #for(ss in 2:length(tissues)){
# study<-tissues[ss]
# print(study)
# 
# 
# wasp_1<- fread(wasp$file_name[wasp$full_name==study][1])
# colnames(wasp_1)[1:2]<- c("chr", "start")
# 
# sam<-unique(wasp_1$SAMPLE_ID)[1]
# wasp_1_sam<-wasp_1 %>% filter(SAMPLE_ID==sam) 
# 
# 
# ase_df<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==sam][1]) %>%
#   filter(pred_genotype==2, coverage>=8) %>% 
#   mutate(ref_ratio=ref_count/coverage)
# 
# 
# ase_wasp<- inner_join(ase_df,wasp_1_sam)
# 
# 
# 
# 
# no_wasp_1<- fread(no_wasp$file_name[no_wasp$full_name==study][1])
# colnames(no_wasp_1)[1:2]<- c("chr", "start")
# 
# no_wasp_1_sam<-wasp_1 %>% filter(SAMPLE_ID==sam) 
# ase_no_wasp<- inner_join(ase_df,no_wasp_1_sam)
# 
# 
# 
# join_ase<- inner_join(wasp_1,no_wasp_1, by=c("chr", "start","SAMPLE_ID"))
# rm(no_wasp_1)
# rm(wasp_1)
# 
# plot_df<- inner_join(ase_wasp,ase_no_wasp,by=c("chr", "start","SAMPLE_ID"))%>% 
#   mutate(alt_count.x=coverage.x-ref_count.x,alt_count.y=coverage.y-ref_count.y,
#          mean_ref=(log2(ref_count.x)+log2(ref_count.y))/2,
#          ratio_ref= log2(ref_count.x)-log2(ref_count.y),
#          mean_alt=(log2(alt_count.x)+log2(alt_count.y))/2,
#          ratio_alt= log2(alt_count.x)-log2(alt_count.y),
#          mean_ref_ratio=(log2(ref_ratio.x)+log2(ref_ratio.y))/2,
#          ratio_ref_ratio= log2(ref_ratio.x)-log2(ref_ratio.y))
# 
# 
# plot_df<-join_ase %>% filter(SAMPLE_ID== unique(join_ase$SAMPLE_ID)[1]) %>% 
#   mutate(mean_ref=(log2(REF_COUNT.x)+log2(REF_COUNT.y))/2,
#          ratio_ref= log2(REF_COUNT.x)-log2(REF_COUNT.y),
#          mean_alt=(log2(ALT_COUNT.x)+log2(ALT_COUNT.y))/2,
#          ratio_alt= log2(ALT_COUNT.x)-log2(ALT_COUNT.y),
#          mean_ref_ratio=(log2(REF_RATIO.x)+log2(REF_RATIO.y))/2,
#          ratio_ref_ratio= log2(REF_RATIO.x)-log2(REF_RATIO.y))
# 
# pdf(file="~/plot/ASE/waspVSnowasp_scatter_recount.pdf", width = 10, height = 6)
# 
# 
# ggplot(plot_df, aes(x=mean_ref, y=ratio_ref))+
#   geom_point(alpha=0.4)+
#   geom_abline(slope=0,color="red")+
#   labs(y="log2 ratio wasp_ref over no_wasp_ref",
#        x= "mean ref counts",
#        title="MA plot of ref counts for one sample")
# 
# ggplot(plot_df, aes(x=mean_alt, y=ratio_alt))+
#   geom_point(alpha=0.4)+
#   geom_abline(slope=0,color="red")+
#   labs(y="log2 ratio wasp_alt over no_wasp_alt",
#        x= "mean alt counts",
#        title="MA plot of alt counts for one sample")
# 
# 
# ggplot(plot_df, aes(x=REF_RATIO.x, y=REF_RATIO.y))+
#   geom_point(alpha=0.4)+
#   geom_abline(slope=1,color="red")+
#   labs(y="wasp ref_ratio",
#        x= "no-wasp ref_ratio",
#        title="ref_ratio comparison for one sample")
# dev.off()
# 

#======================================================

no_wasp_1_1<-no_wasp_1 %>% 
  filter(LOW_MAPABILITY<1,MAPPING_BIAS_SIM<1,GENOTYPE_WARNING<1) %>% 
  select(chr,start,BINOM_P_ADJUSTED,SAMPLE_ID)
wasp_1_1<-wasp_1 %>% 
  filter(LOW_MAPABILITY<1,MAPPING_BIAS_SIM<1,GENOTYPE_WARNING<1)%>% 
  select(chr,start,BINOM_P_ADJUSTED,SAMPLE_ID)


join_ase<- inner_join(wasp_1_1,no_wasp_1_1, by=c("chr", "start","SAMPLE_ID"))

join_ase<-join_ase %>% 
  group_by(SAMPLE_ID) %>% 
  summarize(sig_snps_wasp=sum(BINOM_P_ADJUSTED.x<0.05,na.rm=T),
            sig_snps_no_wasp=sum(BINOM_P_ADJUSTED.y<0.05,na.rm=T),
            sig_snps_both=sum(BINOM_P_ADJUSTED.x<0.05 & BINOM_P_ADJUSTED.y<0.05 ,na.rm=T))

join_ase$diff_wap_nowasp<-join_ase$sig_snps_no_wasp-join_ase$sig_snps_wasp






