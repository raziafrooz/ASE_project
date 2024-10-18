setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(ggplot2)

#########################
#Plot data analysis: Liver
#########################
tissue<-"Liver"

load(paste0("~/hansen_lab/ASE/test_ASE/", tissue, "_ase.rda") )
true_gtex<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/true_gtex_", tissue, ".rds"))
colnames(true_gtex)[1:2]<-c("chr","pos")

sam<- unique(ase_all$sample_id)[2]
ase_all_1 <- ase_all %>% filter(sample_id %in% sam)
true_gtex_1 <- true_gtex %>% filter(sample_id %in% sam)

ase_all_1$aFC<-log2((ase_all_1$alt + 1)/(ase_all_1$ref + 1))

ase_same<-inner_join(true_gtex_1,ase_all_1, by=c("chr", "pos", "sample_id"))


ase_all_1 %>% filter(q_val<0.01)


#----------------
#Get gene annotations
#------------------
library(AnnotationHub)
ah <- AnnotationHub()
ah <- query(ah, c("v26","GENCODE","Homo sapiens","GRch38")) #v26 was used in GTExV8
TxDb<-ah[["AH75155"]]
gene_grch38<-genes(TxDb,columns=c("TXID", "TXNAME"))
tx_38<-transcripts(TxDb,columns=c("TXNAME","GENEID"))

ase_gr<-makeGRangesFromDataFrame(ase_all_1, seqnames.field="chr",
                                 start.field="pos", end.field="pos")
ov<-findOverlaps(tx_38,ase_gr)

ase_all_1$gene_id<-NA
ase_all_1$gene_id[subjectHits(ov)]<-tx_38$GENEID[queryHits(ov)]
ase_all_1$tx_id<-NA
ase_all_1$tx_id[subjectHits(ov)]<-tx_38$TXNAME[queryHits(ov)]


x<-ase_all_1 %>% group_by(gene_id, sample_id) %>%
  summarize(total=n(), ae_hit= sum(q_val<0.01), aFC_max=max(aFC)) %>% ungroup()

x2<-ase_all_1 %>% group_by(tx_id, sample_id) %>%
  summarize(total=n(), ae_hit= sum(q_val<0.01), aFC_max=max(aFC), gene_id=as.character(gene_id[1])) %>% ungroup()

y2<-y[y$gene_id %in% y$gene_id[duplicated(y$gene_id)][2],]
as.data.frame(y2)[1:6,]
y<-x2 %>% filter(ae_hit>0) %>%  group_by(gene_id) %>% summarize(samp=length(sample_id)/9)

x[x$ae_hit>0,][10:20,]
summary(x[x$ae_hit>0,])
y[y$samp>0.5,]
x %>% filter(ae_hit>0) %>%  group_by(gene_id) %>% summarize(total=n())
pdf(file="~/plot/ASE/Gene_analysis.pdf", width = 10, height = 4)
ggplot(y)+
  geom_histogram(aes(samp))
ggplot(x[x$ae_hit>0,][1:100,])+
  geom_boxplot(aes(x=gene_id, y=aFC_max))
dev.off()

x$gene_id<-unlist(lapply(strsplit(x$gene_id, "[.]"), function(x){
  x[[1]]
}))
#--------------------
shet<-read.csv("~/naturalSelection/tf_selection/data/regeneron.txt", sep="\t")
x$shet<- shet$mean[match(x$gene_id, shet$GeneId)]
x2<-x[x$ae_hit>0,]

Q<-quantile(x2$shet, probs = seq(.1, 1, by = .1),na.rm=T)
x2<-x2 %>% group_by(gene_id) %>% mutate(sample_number=length(sample_id),
                                    shet_decile=case_when(shet <= Q[[1]] ~ "1",
                                                          Q[[1]]< shet & shet<= Q[[2]] ~ "2",
                                                          Q[[2]]< shet & shet<= Q[[3]] ~ "3",
                                                          Q[[3]]< shet & shet<= Q[[4]] ~ "4",
                                                          Q[[4]]< shet & shet<= Q[[5]] ~ "5",
                                                          Q[[5]]< shet & shet<= Q[[6]] ~ "6",
                                                          Q[[6]]< shet & shet<= Q[[7]] ~ "7",
                                                          Q[[7]]< shet & shet<= Q[[8]] ~ "8",
                                                          Q[[8]]< shet & shet<= Q[[9]] ~ "9",
                                                          shet > Q[[9]] ~ "10"))

x_2<-x2[!duplicated(x2$gene_id),]


pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
ggplot(x_2)+
  geom_boxplot(aes(x=shet_decile,y=sample_number))
dev.off()

x_2$perc<-x_2$sample_number/47
x_2[x_2$perc>0.5,]
