setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(rtracklayer)

qc_df<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/gtex_qc_metadata.csv.gz")
qc_df<-qc_df %>% mutate(overlap= (star.average_input_read_length)-bc_frag.mode_length)


path_id<-"/dcl01/hansen/data/arazi/ASE/dbGap/phe000039.v1.GTEx_v8_ASE_WASP.expression-matrixfmt-ase.c1/GTEx_Analysis_v8_ASE_WASP_counts_by_subject/"
tissues_names<-as.data.frame(list.files(path =path_id))
colnames(tissues_names)<-"file_name"
tissues_names$indv<-gsub("\\..*","",tissues_names$file_name)
tissues_names$file_name<-paste0(path_id, tissues_names$file_name)

gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)

#plotFile<-fread("/dcs07/hansen/data/recount_ASE/data/gtexVSrecount_wasp.csv.gz")
#colnames(plotFile)[2]<-"sample_id_rep"
tt<-fread("/dcs07/hansen/data/recount_ASE/data/gtex_recountPipeline.csv.gz")
tt$sample_id_rep<-str_sub(tt$sample_id, end= -3)
#xx<-left_join(tt,plotFile)

tt$overlap<-qc_df$overlap[match(tt$sample_id ,qc_df$external_id)]
tt$overlap_group<-cut_number(tt$overlap,6)
tt$indv<-sapply(strsplit(tt$sample_id_rep,"-"), function(xx){ paste0(xx[1],"-",xx[2])  })

i=which(tt$overlap< -80)[1]
#all<-c()
for(i in 1747:nrow(tt)){
print(i)
indv<-tt$indv[i]
rep<-tt$sample_id_rep[i]
ov<-tt$overlap[i]
if(length(tissues_names$file_name[tissues_names$indv==indv])>0){
gtex_ase<- fread(tissues_names$file_name[tissues_names$indv==indv])
colnames(gtex_ase)[1:2]<- c("chr", "start")
gtex_ase_1<-gtex_ase %>% filter(SAMPLE_ID ==rep)# %>% select(c(chr, start,REF_COUNT,ALT_COUNT,TOTAL_COUNT,REF_RATIO))

recount_ase<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==rep][1]) %>%
  filter(pred_genotype==2, coverage>=8) %>% 
  mutate(alt_count=coverage-ref_count,ref_ratio=ref_count/coverage) %>% 
  select(-c(pred_genotype,predicted_accuracy, M,S,AF))



joined<-inner_join(gtex_ase_1,recount_ase) %>% 
  mutate(ratio_total=log2(TOTAL_COUNT) - log2(coverage),
            mean_total=(log2(TOTAL_COUNT) + log2(coverage))/2,
            ratio_ref=log2(REF_COUNT) - log2(ref_count),
            mean_ref=(log2(REF_COUNT) + log2(ref_count))/2,
            ratio_alt=log2(ALT_COUNT) - log2(alt_count),
            mean_alt=(log2(ALT_COUNT) + log2(alt_count))/2)


one<-data.frame(overlap=ov,
           ref_diff=median(joined$ratio_ref),
           alt_diff=median(joined$ratio_alt),
           total_diff=median(joined$ratio_total))
all<-rbind(all, one)
}
}
#fwrite(all, "~/test/gtex_qc.csv.gz")
joined<-joined[-which(joined$LOW_MAPABILITY>0 | joined$MAPPING_BIAS_SIM>0 | joined$GENOTYPE_WARNING>0),]

pdf(file="~/plot/ASE/indv_no_ov.pdf", width = 10, height = 6)

pp=ggplot(joined, aes(y=ratio_total, x=mean_total))+
  geom_point(alpha=0.3)+
  geom_hline(yintercept = 0,linetype = "dashed", alpha=0.5, color="red")+
  geom_smooth()+
  labs(title=paste0("sample has overlap = ", ov),
       subtitle=rep)
print(pp)

pp=ggplot(joined, aes(y=ratio_ref, x=mean_ref))+
  geom_point(alpha=0.3)+
  geom_hline(yintercept = 0,linetype = "dashed", alpha=0.5, color="red")+
  geom_smooth()+
  labs(title=paste0("sample has overlap = ", ov),
       subtitle=rep)

print(pp)

pp=ggplot(joined, aes(y=ratio_alt, x=mean_alt))+
  geom_point(alpha=0.3)+
  geom_hline(yintercept = 0,linetype = "dashed", alpha=0.5, color="red")+
  geom_smooth()+
  labs(title=paste0("sample has overlap = ", ov),
       subtitle=rep)

print(pp)
dev.off()

all$overlap_group<-cut_number(all$overlap,5)

sp<-pivot_longer(all, !c(overlap,overlap_group), names_to = "difference", values_to = "value")

pdf(file="~/plot/ASE/count_summary.pdf", width = 10, height = 6)

pp=ggplot(sp, aes(x=difference, y=value, color=overlap_group))+
  geom_boxplot()+
  labs(title="Summary plot:Median of ratios",
       subtitle="negative values indicate that Recount has higher counts than GTEx",
       y="Median of ratio",
       x="Count group")

print(pp)
dev.off()
