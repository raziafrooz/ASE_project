library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)

qc_df<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/gtex_qc_metadata.csv.gz")
qc_df<-qc_df %>% mutate(overlap= (star.average_input_read_length)-bc_frag.mode_length)
gtex_uni_norm<-fread("/dcs07/hansen/data/recount_ASE/data/gtex_uniNorm.csv")
colnames(gtex_uni_norm)[2]<-"SAMPLE_ID"
qc_df$SAMPLE_ID<-str_sub(qc_df$external_id, end= -3)

#============================================================
#Read in wasp indv:
path_id<-"/dcl01/hansen/data/arazi/ASE/dbGap/phe000039.v1.GTEx_v8_ASE_WASP.expression-matrixfmt-ase.c1/GTEx_Analysis_v8_ASE_WASP_counts_by_subject/"
tissues_names<-as.data.frame(list.files(path =path_id))
colnames(tissues_names)<-"file_name"
tissues_names$indv<-gsub("\\..*","",tissues_names$file_name)
tissues_names$file_name<-paste0(path_id, tissues_names$file_name)
#============================================================



gtex_w_ov<-qc_df %>% filter(SMGEBTCHT=="TruSeq.v1",overlap>40) %>% 
  left_join(gtex_uni_norm,by=c("SAMPLE_ID")) %>% 
  filter(uni_norm_mean<0.08)
row_id=2
gtex_w_ov$adj_ref_ratio<-NA
adj_stat_all<-c()
for(row_id in 1:nrow(gtex_w_ov)){
sample_id<-gtex_w_ov$SAMPLE_ID[row_id]
study<-gtex_w_ov$study[row_id]
ov<-gtex_w_ov$overlap[row_id]
indv<-gtex_w_ov$SUBJID[row_id]
print(row_id)

ase_df<-fread(gtex_metadata$genotypedSamples[which(gtex_metadata$sample_id_rep==sample_id)][1]) %>% 
  filter(pred_genotype==2, coverage>=8) %>% 
  mutate(ref_ratio=ref_count/coverage,
         alt_count=coverage-ref_count,
         ratio=log2(ref_count/alt_count),
         mean=(log2(ref_count)+log2(alt_count))/2)

#Fix the counts by adjusting the MA plot to mean around 0
ratio_adj<-median(ase_df$ratio)
if(ratio_adj!=0){
  
  adj<-ase_df$ratio-ratio_adj
  ase_df$adj_alt<-round(ase_df$coverage/((2^adj)+1),0)

stopifnot(sum(is.na(ase_df$adj_alt))==0)

ase_df$adj_ref<-ase_df$coverage-ase_df$adj_alt
print(median(ase_df$adj_ref/ase_df$coverage))
gtex_w_ov$adj_ref_ratio[row_id]<-median(ase_df$adj_ref/ase_df$coverage)
#------------------------
#compare with counts in WASP:
wasp<-tissues_names$file_name[tissues_names$indv==indv]
if(length(wasp)>0){
gtex_ase<- fread(wasp)
colnames(gtex_ase)[1:2]<- c("chr", "start")
gtex_ase_1<-gtex_ase %>% filter(SAMPLE_ID ==sample_id) %>% select(c(chr, start,REF_COUNT,ALT_COUNT,TOTAL_COUNT,REF_RATIO))

joined<-ase_df %>% 
  dplyr::select(c(chr, start,adj_alt, adj_ref,coverage,ref_count,alt_count)) %>% 
  inner_join(gtex_ase_1,by = c("chr", "start")) 


adj_stat<-data.frame(overlap=ov,
                     sample_id,
                     ref_before_adj=median((joined$ref_count -joined$REF_COUNT)/joined$ref_count),
           ref_after_adj=median((joined$adj_ref -joined$REF_COUNT)/joined$adj_ref),
           alt_before_adj=median((joined$alt_count -joined$ALT_COUNT)/joined$alt_count),
           alt_after_adj=median((joined$adj_alt -joined$ALT_COUNT)/joined$adj_alt),
           total_diff=median((joined$coverage -joined$TOTAL_COUNT)/joined$coverage))


adj_stat_all<-rbind(adj_stat_all,adj_stat)


}
}
}

#fwrite(adj_stat_all, file="~/test/gtex_fixed.csv")
adj_stat_all<-fread("~/test/gtex_fixed.csv")

plot_df<-adj_stat_all %>% dplyr::select(overlap,sample_id,alt_before_adj,alt_after_adj) %>% pivot_longer(.,
                      !c("overlap","sample_id"), names_to="allele", values_to="count_difference")
plot_df$overlap_group<-cut_number(plot_df$overlap, 5)
pdf(file="~/plot/ASE/adjustVSwasp.pdf", width = 10, height = 6)
ggplot(data=plot_df,aes(x=overlap_group, y=count_difference, fill=allele))+
  geom_boxplot()+
  scale_fill_brewer(palette="BuPu")+
  labs(title="After correcting overlap issue in 4,367 GTEx samples")
dev.off()

