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

#-----------------------------------

seq_mean=seq(0,15,by=0.01)


gtex_w_ov<-qc_df %>% filter(SMGEBTCHT=="TruSeq.v1", overlap>0) %>% 
  left_join(gtex_uni_norm,by=c("SAMPLE_ID")) %>% 
  filter(uni_norm_mean<0.08)

gtex_w_ov$adj_ref_ratio<-NA
for(row_id in 1:nrow(gtex_w_ov)){
sample_id<-gtex_w_ov$SAMPLE_ID[row_id]
study<-gtex_w_ov$study[row_id]
ov<-gtex_w_ov$overlap[row_id]
print(row_id)

ase_df<-fread(gtex_metadata$genotypedSamples[which(gtex_metadata$sample_id_rep==sample_id)][1]) %>% 
  filter(pred_genotype==2, coverage>=8) %>% 
  mutate(ref_ratio=ref_count/coverage,
         alt_count=coverage-ref_count,
         ratio=log2(ref_count/alt_count),
         mean=(log2(ref_count)+log2(alt_count))/2)



ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)


zz<-ase_df %>% dplyr::select(chr,start,coverage, ref_count,alt_count, ratio,mean, mean_20tile) %>% 
  group_by(mean_20tile) %>% reframe(median_ratio=median(ratio),
                                    mean_ratio=mean(ratio)) %>% 
  ungroup() %>% 
  mutate(ratio_adj=round(mean_ratio- 0,3))


id<-which(zz$ratio_adj!=0)
ase_df$adj_alt<-NA
for(kk in 1:length(id)){

  i<-id[kk]
  group<-zz$mean_20tile[i]
  ratio_adj<-zz$ratio_adj[i]
  
  indx<-which(ase_df$mean_20tile==group)
  adj<-ase_df$ratio[indx]-ratio_adj
  ase_df$adj_alt[indx]<-round(ase_df$coverage[indx]/((2^adj)+1),0)
}

ase_df$adj_alt[is.na(ase_df$adj_alt)]<-ase_df$alt_count[is.na(ase_df$adj_alt)]

#ase_df[,-c("mean_20tile","gtex_median",  "gtex_mean","chr","start")]



ase_df$adj_ref<-ase_df$coverage-ase_df$adj_alt

print(median(ase_df$adj_ref/ase_df$coverage))
gtex_w_ov$adj_ref_ratio[row_id]<-median(ase_df$adj_ref/ase_df$coverage)
}







pdf(file="~/plot/ASE/MA_test2.pdf", width = 10, height = 6)
ggplot(data=ase_df,aes(y=log2(ref_count), x=log2(alt_count)))+
  geom_point(alpha=0.4)+
  geom_abline(slope=1, color="red")+
  labs(title=paste0("Before adjustment:", sample_id),
       subtitle=paste0("overlap is = ", ov))


ggplot(data=ase_df,aes(y=log2(adj_ref), x=log2(adj_alt)))+
  geom_point(alpha=0.5)+
  geom_abline(slope=1, color="red")+
  labs(title=paste0("After adjustment:", sample_id),
       subtitle=paste0("overlap is = ", ov))


ggplot(data=ase_df,aes(y=log2(ref_count), x=log2(alt_count)))+
  geom_point(alpha=0.5)+
  geom_abline(slope=1, color="red")+
  geom_point(data=ase_df,aes(y=log2(adj_ref), x=log2(adj_alt)),color="blue",alpha=0.5)


dev.off()




# for(k in 1:nrow(qc_df_good)){
#   sample_id<-qc_df_good$SAMPLE_ID[k]
#   study<-qc_df_good$study[k]
#   print(k)
#   
#     ase_df<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==sample_id][1]) %>% 
#       filter(pred_genotype==2, coverage>=8) %>% 
#       mutate(ref_ratio=ref_count/coverage,
#              alt_count=coverage-ref_count,
#              ratio=log2(ref_count/alt_count),
#              mean=(log2(ref_count)+log2(alt_count))/2)
#     
# 
#     
#     ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)
#     
#     
#     
#     q_line<-ase_df %>% group_by(mean_20tile) %>%
#       reframe(median_ratio=median(ratio),
#               mean_ratio=mean(ratio))
#     q_line$sample<-sample_id
#     q_line$study<-study
#     
#     
#     q_line_geu<-rbind(q_line_geu,q_line)
#     
#   }
#   
# 
# test_line<-q_line_geu %>%
#   group_by(mean_20tile) %>%
#   reframe(gtex_median=median(median_ratio),
#           gtex_mean=median(mean_ratio))
#fwrite(test_line, "~/test/gtex_adj_count.csv")
#test_line<-fread("~/test/gtex_adj_count.csv")

# ase_df<-ase_df %>%
#   dplyr::select(chr,start,coverage, ref_count,alt_count, ratio,mean, mean_20tile) %>%
#   left_join(test_line,by="mean_20tile")


# pdf(file="~/plot/ASE/MA_test.pdf", width = 10, height = 6)
# ggplot(data=ase_df,aes(y=ratio, x=mean))+
#   geom_point(alpha=0.5)+
#   geom_hline(yintercept=0, color="red")+
#   geom_point(data=ase_df,aes(x=mean, y=gtex_mean),color="blue")
# 
# ggplot(data=zz,aes(y=ratio, x=mean))+
#   geom_point(alpha=0.5)+
#   geom_hline(yintercept=0, color="red")+
#   geom_point(data=ase_df,aes(x=mean, y=gtex_mean),color="blue")
# 
# 
# dev.off()

# 
# ase_df %>% filter(mean_20tile== "(1.4,1.41] ") %>% head(8)
# zz<-ase_df %>%
#   group_by(mean_20tile) %>% reframe(median_ratio=median(ratio),
#                                     mean_ratio=mean(ratio)) %>%
#   left_join(test_line,by="mean_20tile") %>% ungroup() %>%
#   mutate(ratio_adj=round(mean_ratio-gtex_mean,3))





median(ase_df$ref_count/ase_df$coverage)

median(ase_df$adj_ref/ase_df$coverage)
ase_df[which(log2(ase_df$adj_alt)<3),][1:4,]
