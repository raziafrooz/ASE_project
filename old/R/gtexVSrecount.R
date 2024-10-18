setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)
tissue<-"Lung"

load(paste0("~/hansen_lab/ASE/test_ASE/", tissue, "_ase.rda") ) #ase_all
true_gtex<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/true_gtex_", tissue, ".rds"))
colnames(true_gtex)[1:2]<-c("chr","pos")


sum_G<-true_gtex %>% group_by(sample_id) %>%  summarize(m_ratio_G=median(REF_RATIO)) %>%  arrange(desc(m_ratio_G))

sum_R<-ase_all %>% group_by(sample_id) %>%  summarize(m_ratio_R=median(ref_ratio)) %>%  arrange(desc(m_ratio_R))

tail(sum_G)
tail(sum_R)

#in recount we have this sample with low ref_ratio GTEX-NPJ8-0326-SM-2D7VV.1 = 0.458 vs GTEx GTEX-NPJ8-0326-SM-2D7VV.1 = 0.5
# only focus on this sample for now
library(cowplot)

sam<- sum_R$sample_id
for (i in 1:nrow(sum_R)){
  sam<- sum_R$sample_id[i]
ase_all_1 <- ase_all %>% filter(sample_id %in% sam)
true_gtex_1 <- true_gtex %>% filter(sample_id %in% sam)

ase_all_2<-ase_all_1 %>% left_join(true_gtex_1) %>%  mutate(in_gtex= case_when(!is.na(TOTAL_COUNT) ~ "yes" ,is.na(TOTAL_COUNT) ~ "no"))

#sum(ase_all_2$in_gtex =="yes") #2033 from Recount3 is in GTEx
##sum(ase_all_2$in_gtex =="no") #184 from Reocunt3 are not in GTEx
#sum(ase_all_2$in_gtex =="no" & ase_all_2$true_genotype!=2) #Out of 184 dicrepencies 99 are genotyping error in Reocunt3

#ase_all_2[ase_all_2$in_gtex=="yes",c("total","ref", "alt","TOTAL_COUNT", "REF_COUNT", "ALT_COUNT")]
ase_all_2<- ase_all_2 %>% mutate(ratio_ref=log2(REF_COUNT) - log2(ref),
                     mean_ref=(log2(REF_COUNT) + log2(ref))/2,
                     ratio_alt=log2(ALT_COUNT) - log2(alt),
                     mean_alt=(log2(ALT_COUNT) + log2(alt))/2,
                     ratio_total=log2(TOTAL_COUNT) - log2(total),
                     mean_total=(log2(TOTAL_COUNT) + log2(total))/2)


pdf(file=paste0("~/plot/ASE/GTEx_vs_recount/",sam,".pdf"), width = 10, height = 6)

p= ggplot(ase_all_2[ase_all_2$in_gtex=="yes",]) + 
  geom_histogram(aes(REF_RATIO), fill="darkblue", alpha=0.3)+
  geom_histogram(aes(ref_ratio), fill="salmon", alpha=0.3)+
  labs(title="Same SNPs that are in Recount3: blue is GTEx, pink is Recount")+
   annotate("text", x= 0.7, y=230, label= paste0("total gtex=" , nrow(true_gtex_1), ", total recount=", nrow(ase_all_1)))+
   annotate("text", x= 0.7, y=220, label= paste0("gtex ref-ratio=" , round(median(ase_all_2$REF_RATIO,na.rm=T),3), ", recount ref-ratio=",round(median(ase_all_2$ref_ratio,na.rm=T),3)))+
   annotate("text", x= 0.7, y=200, label= paste0(sum(ase_all_2$in_gtex =="yes")," from Recount3 is in GTEx")) + 
   annotate("text", x = 0.7, y=180, label = paste0(sum(ase_all_2$in_gtex =="no"), " from Recount3 is not is GTEx"))+
   annotate("text", x = 0.7, y=160, label = paste0("out of ", sum(ase_all_2$in_gtex =="no"), 
                                                  " discrepancies ", sum(ase_all_2$in_gtex =="no" & ase_all_2$true_genotype!=2), " are genotyping error in Recount3"))
print(p)
p= ggplot(ase_all_2[ase_all_2$in_gtex=="yes",])+
   geom_point(aes(x=log10(ALT_COUNT), y=log10(REF_COUNT)), color="darkblue", alpha=0.3)+
   geom_point(aes(x=log10(alt), y=log10(ref)), color="salmon", alpha=0.3)
print(p)

 p1 <- ggplot(ase_all_2[ase_all_2$in_gtex=="yes",]) + 
   geom_point(aes(x=log10(ALT_COUNT), y=log10(REF_COUNT)), color="lightblue")
 p2 <- ggplot(ase_all_2[ase_all_2$in_gtex=="yes",]) + 
   geom_point(aes(x=log10(alt), y=log10(ref)), color="salmon")
 
 p=plot_grid(p1, p2, labels = c('GTEx', 'Recount'))
 print(p)
 

p1 <- ggplot(ase_all_2[ase_all_2$in_gtex=="yes",]) + 
  geom_boxplot(aes(x=sample_id, y=REF_RATIO), color="lightblue")
p2 <- ggplot(ase_all_2[ase_all_2$in_gtex=="yes",]) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio), color="salmon")

p=plot_grid(p1, p2, labels = c('GTEx', 'Recount'))
print(p)

p=ggplot(ase_all_2[ase_all_2$in_gtex=="yes",],aes(y=ratio_total,x=mean_total))+
  geom_point(alpha=0.4)+
  geom_smooth()
print(p)

p=ggplot(ase_all_2[ase_all_2$in_gtex=="yes",],aes(y=ratio_alt,x=mean_alt))+
  geom_point(alpha=0.4)+
  geom_smooth()
print(p)

p=ggplot(ase_all_2[ase_all_2$in_gtex=="yes",],aes(y=ratio_ref,x=mean_ref))+
  geom_point(alpha=0.4)+
  geom_smooth()
print(p)

p=ggplot(ase_all_2[ase_all_2$in_gtex=="yes",])+
  geom_point(aes(x=REF_RATIO,y=ref_ratio), alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="salmon")
print(p)

dev.off()
}

