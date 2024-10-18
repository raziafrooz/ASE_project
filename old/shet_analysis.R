setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(ggplot2)


#----------------
#Get gene annotations
#------------------
library(AnnotationHub)
ah <- AnnotationHub()
ah <- query(ah, c("v26","GENCODE","Homo sapiens","GRch38")) #v26 was used in GTExV8
TxDb<-ah[["AH75155"]]
gene_grch38<-genes(TxDb,columns=c("TXID", "TXNAME"))
tx_38<-transcripts(TxDb,columns=c("TXNAME","GENEID"))


#########################
#Plot data analysis: Liver
#########################


tissues<-c("Liver", "Lung", "Stomach", "Pancreas")

for(k in 1:length(tissues)){
  print(tissues[k])
  tissue<-tissues[k]
load(paste0("~/hansen_lab/ASE/test_ASE/", tissue, "_ase.rda") )

#sam<- unique(ase_all$sample_id)[2]
#ase_all_1 <- ase_all %>% filter(sample_id %in% sam)
ase_all<-ase_all %>% group_by(sample_id) %>% filter(!(median(ref_ratio) < 0.4), !(median(ref_ratio) > 0.6) )
ase_all$aFC<-log2((ase_all$alt + 1)/(ase_all$ref + 1))

#----------------
#Get gene annotations
#------------------
ase_gr<-makeGRangesFromDataFrame(ase_all, seqnames.field="chr",
                                 start.field="pos", end.field="pos")
ov<-findOverlaps(tx_38,ase_gr)

ase_all$gene_id<-NA
ase_all$gene_id[subjectHits(ov)]<-as.character(tx_38$GENEID[queryHits(ov)])

print("Selecting SNP-wise ASE")
ase_all<-ase_all %>%  group_by(gene_id,sample_id) %>% filter(total==max(total)) %>% ungroup()

if(k==1){
  x<-ase_all %>% group_by(gene_id) %>%
    summarize(total=n(), ae_hit= sum(q_val<0.05)) %>% ungroup()
}else{
  
  x2<-ase_all %>% group_by(gene_id) %>%
    summarize(total=n(), ae_hit= sum(q_val<0.05)) %>% ungroup()
  id<-match(x$gene_id,x2$gene_id)
  tt<-x2$total[id]
  tt[is.na(tt)]<-0
  
  ase<-x2$ae_hit[id]
  ase[is.na(ase)]<-0
  
  x$total<-x$total+tt
  x$ae_hit<-x$ae_hit+ase
}
}
xx<-x



xx$percent<-xx$ae_hit/xx$total
# summary(x[x$ae_hit>0,])
# y[y$samp>0.5,]
# x %>% filter(ae_hit>0) %>%  group_by(gene_id) %>% summarize(total=n())
# pdf(file="~/plot/ASE/Gene_analysis.pdf", width = 10, height = 4)
# ggplot(y)+
#   geom_histogram(aes(samp))
# ggplot(x[x$ae_hit>0,][1:100,])+
#   geom_boxplot(aes(x=gene_id, y=aFC_max))
# dev.off()
# 
 xx$gene_id<-unlist(lapply(strsplit(xx$gene_id, "[.]"), function(y){
   y[[1]]
 }))
#--------------------
shet<-read.csv("~/naturalSelection/tf_selection/data/regeneron.txt", sep="\t")
xx$shet<- shet$mean[match(xx$gene_id, shet$GeneId)]


Q<-quantile(xx$shet, probs = seq(.1, 1, by = .1),na.rm=T)
Q<-quantile(xx$shet,na.rm=T)
xx2<-xx %>% mutate(shet_decile=case_when(shet <= Q[[2]] ~ "1",
  Q[[2]]< shet & shet<= Q[[3]] ~ "2",
  Q[[3]]< shet & shet<= Q[[4]] ~ "3",
  Q[[4]]< shet ~ "4"))

xx2<-xx %>% mutate(shet_decile=case_when(#shet <= Q[[1]] ~ "1",
                   Q[[1]]< shet & shet<= Q[[2]] ~ "2",
                   Q[[2]]< shet & shet<= Q[[3]] ~ "3",
                   Q[[3]]< shet & shet<= Q[[4]] ~ "4",
                   Q[[4]]< shet ~ "5"))
                   
                   # Q[[5]]< shet & shet<= Q[[6]] ~ "6",
                   # Q[[6]]< shet & shet<= Q[[7]] ~ "7",
                   # Q[[7]]< shet & shet<= Q[[8]] ~ "8",
                   # Q[[8]]< shet & shet<= Q[[9]] ~ "9",
                   # shet > Q[[9]] ~ "10"))

x_2<-x2[!duplicated(x2$gene_id),]

library(ggridges)

pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
ggplot(xx2 %>% filter(!is.na(shet_decile)))+
  geom_density(aes(percent, color=shet_decile,fill=shet_decile), alpha=0.6)+lims(x=c(0,0.3))
ggplot(xx2 %>% filter(!is.na(shet_decile)), aes(x = percent, y = shet_decile)) + geom_density_ridges2()
dev.off()

x_2$perc<-x_2$sample_number/47
x_2[x_2$perc>0.5,]
