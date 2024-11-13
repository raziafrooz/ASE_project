library(AnnotationHub)
library("org.Hs.eg.db")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

#===========================================================
tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("~/ASE-data/data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

tissues_names$file_name<-paste0("/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/", tissues_names$file_name)

gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)

#-----------------------------------------------------------------
#-----------------------------------------------------------------
xx<-read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/metadata/all_gtex_metadata.csv")
gtex_sim<-fread("/dcs07/hansen/data/recount_ASE/data/gtex_simulation.csv.gz")
gtex_sim_gr<-makeGRangesFromDataFrame(gtex_sim,seqnames="chr",start.field ="start",end.field = "start")


qc_df<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/gtex_qc_metadata.csv.gz")
qc_df<-qc_df %>% mutate(overlap= (star.average_input_read_length)-bc_frag.mode_length)
qc_df$SAMPLE_ID<-str_sub(qc_df$external_id, end= -3)
qc_df<-qc_df %>% 
  mutate(overlap= (star.average_input_read_length)-bc_frag.mode_length) %>% 
  filter(SMGEBTCHT=="TruSeq.v1", overlap<200)%>% 
  dplyr::select(SAMPLE_ID,overlap )
# 
# nrow(xx)
# for (k in 2:3){
#   print(k)
#   study<-unique(xx$study)[k]
#   print(study)
#   wasp_1<- fread(tissues_names$file_name[tissues_names$full_name==study][1]) %>% 
#     filter(LOW_MAPABILITY<1,MAPPING_BIAS_SIM<1,GENOTYPE_WARNING<1)
#   colnames(wasp_1)[1:2]<- c("chr", "start")
#   
# 
#   sum_one<-wasp_1 %>% group_by(SAMPLE_ID) %>% summarize(total_sum=sum(TOTAL_COUNT),
#                                                         ref_sum=sum(REF_COUNT),
#                                                         Nsnp=n())
#   sum_one$study<-study
#   
#   sum_total<-rbind(sum_total,sum_one)
# }
# #back<-sum_total
# try1<-right_join(qc_df,sum_total)
# 
# try1<-right_join(try1,recount_sig_snps)
# 
# pdf(file="~/plot/ASE/test2.pdf", width = 10, height = 6)
# 
#   
#   pp=ggplot(try1,aes(y=total_sum, x=overlap))+
#     geom_point()
#   
#   print(pp)
#   
#   pp=ggplot(try1,aes(y=total_sum, x=sig_gtex))+
#     geom_point()
#   
#   print(pp)
#   
#   pp=ggplot(try1,aes(y=Nsnp, x=overlap))+
#     geom_point()
#   
#   print(pp)
# 
#   pp=ggplot(try1,aes(y=Nsnp, x=sig_gtex))+
#     geom_point()
#   
#   print(pp)
#   
# dev.off()
# 

#--------------------------------------
#joined_all<-c()

for (k in 9:length(unique(xx$study))){
  study<-unique(xx$study)[k]
  print(study)
  wasp_1<- fread(tissues_names$file_name[which(tissues_names$full_name==study)]) %>% 
    filter(LOW_MAPABILITY<1,MAPPING_BIAS_SIM<1,GENOTYPE_WARNING<1)
  colnames(wasp_1)[1:2]<- c("chr", "start")
  
  
  xx_one<-xx[xx$study==study,]
  
  #for(ss in 3:length(unique(xx$sample_id))){
  for(ss in 1:length(unique(wasp_1$SAMPLE_ID))){
    print(ss)

sam<-unique(wasp_1$SAMPLE_ID)[ss]
sam_id<-xx_one$sample_id_rep[xx_one$sample_id==sam][1]


wasp_1_sam<-wasp_1 %>% filter(SAMPLE_ID==sam) %>% select("chr","start","TOTAL_COUNT","ALT_COUNT","REF_COUNT")
if(nrow(wasp_1_sam)>1){
  
  
  
  ase_df<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==sam][1]) %>%
    filter(pred_genotype==2, coverage>=8) %>% 
    mutate(ref_ratio=ref_count/coverage,
           alt_count=coverage-ref_count) %>% select("chr","start","coverage","alt_count","ref_count")
  wasp_1_sam<-wasp_1_sam %>% select("chr","start","TOTAL_COUNT","ALT_COUNT","REF_COUNT")
  
  
jj<-inner_join(wasp_1_sam,ase_df)%>%
  reframe(total_diff=((coverage-TOTAL_COUNT)/coverage),
          alt_diff=((alt_count-ALT_COUNT)/alt_count))



jj<-jj %>% reframe(quantile_group=c("0","0.25","0.5","0.75","1"),
               q_total=quantile(total_diff),
               q_alt=quantile(alt_diff),
               SAMPLE_ID=sam,
               study)

joined_all<-rbind(joined_all,jj)


}
}
}

#fwrite(joined_all, "~/test/joined_all.csv.gz")
try1<-joined_all
try1$overlap<-qc_df$overlap[match(try1$SAMPLE_ID,qc_df$SAMPLE_ID)]
try1$overlap_group<-cut_number(try1$overlap,7)


pdf(file="~/plot/ASE/Recount_alt_diff.pdf", width = 10, height = 6)

ggplot(data=try1 %>% filter(quantile_group=="0.5"), aes(x=overlap_group, y=q_alt, color=quantile_group)) +
  geom_boxplot()+
  labs(title="on average, what's the difference in alt counts",
       subtitle="recount_alt-gtex_alt/recount_alt",
       y="alt count difference")
  
ggplot(data=try1 %>% filter(quantile_group=="0.5"), aes(x=overlap_group, y=q_total, color=quantile_group)) +
  geom_boxplot()+
  labs(title="on average, what's the difference in total counts",
       subtitle="recount_total-gtex_total/recount_total",
       y="total count difference")
dev.off()  
