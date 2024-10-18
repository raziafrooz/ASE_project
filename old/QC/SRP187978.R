setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)


sra_geno<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")
sra_geno<-sra_geno[sra_geno$study=="SRP187978",]

#Penn cohort2
id<-which(sra_geno$sample_id %in% c("SRR8701271","SRR8701272","SRR8701273",
                                    "SRR8701274","SRR8701275","SRR8701276"))

id<-which(sra_geno$sample_id %in% c("SRR8701277","SRR8701278",
                                    "SRR8701279","SRR8701280",
                                    "SRR8701281","SRR8701282",
                                    "SRR8701283","SRR8701284",
                                    "SRR8701285","SRR8701286"))

#Penn cohort1:
id<-which(sra_geno$sample_id %in% c("SRR8702455","SRR8702456",
                                    "SRR8702457","SRR8702458",
                                    "SRR8702459","SRR8702460",
                                    "SRR8702461","SRR8702462"))


for(k in id){
  print(k)
  sample_id<-sra_geno$sample_id[k]
  study<-sra_geno$study[k]
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda")))
  {
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda"))
    ase_df$sample_id<-sample_id
    ase_df$study<-study
    if(k==id[1]){
      ase_all<-ase_df
    } else{
    ase_all<-rbind(ase_all,ase_df)
    }}}
#ase_all$indv<-3
#ase_all2<-rbind(ase_all2,ase_all)

ase_all %>% group_by(sample_id) %>% 
  summarise(total_snp=n(),q25= quantile(total)[[2]],
            q75= quantile(total)[[4]], avg_cov=mean(total),
            ref_ratio=mean(ref_ratio))

plot<-ase_all2 %>% 
  group_by(chr,pos, indv) %>%
  summarize(ref=sum(ref), alt=sum(alt), total=sum(total)) %>% ungroup()
plot$indv<-as.factor(plot$indv)



pdf(file="~/plot/ASE/SRR8701276_sample3.pdf", width = 10, height = 4)
p=ggplot(ase_all)+
  geom_boxplot(aes(x=sample_id, y=ref_ratio),outlier.shape = NA)+
  geom_boxplot(data= plot, aes(x="combined", y=ref_ratio), fill="salmon",outlier.shape = NA)+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  labs(title="SRR8701276 ref ratio individual runs vs combined (Penn cohort1)")
print(p)

  p=ggplot(plot)+
    geom_point(aes(x=total, y=mine_allele))+xlim(0,100)+ylim(0,100)+
    labs(title="SRR8701276 one sample: mono-allele check (Penn cohort1)")
  print(p)

dev.off()

pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
ggplot()+
  geom_point(data= plot %>% filter(total>=8) %>% summarise(n=n()), aes(x=">=8", y=n))+
  geom_point(data= plot %>% filter(total>=16) %>% summarise(n=n()), aes(x=">=16", y=n))+
  geom_point(data= plot %>% filter(total>=30) %>% summarise(n=n()), aes(x=">=30", y=n))+
  geom_point(data= plot %>% filter(total>=60) %>% summarise(n=n()), aes(x=">=60", y=n))
dev.off()

#-----------------------
#Geuvadis
#-----------------------
geu<-read.csv("~/ASE/data/Geuvadis_metadata.csv")
study<-geu$study[k]
for(t in 1:3){
  print(t)
load(paste0("~/hansen_lab/ASE/test_ASE/", geu$sample_id[t], "_ase.rda"))#named ase_all
if(t==1){
geuvadis<-ase_all
} else {
geuvadis<-rbind(geuvadis,ase_all)
}}
t=1
geuvadis %>% rowwise() %>% mutate(min_allele=min(alt,ref)) %>%  group_by(sample_id) %>% 
  summarise(total_snp=length(total),q25= quantile(total)[[2]], q75= quantile(total)[[4]],
            median_cov=median(total), med_mono=median(min_allele))
plot_geu<-geuvadis[geuvadis$sample_id=="ERR188482",]
median(geuvadis$total)

for(i in c(8,16,30,60)){
  if(i==8){
plot_geu<-geuvadis %>% group_by(sample_id) %>%  filter(total>=i) %>% summarise(n())}else{
  rr<-geuvadis %>% group_by(sample_id) %>%  filter(total>=i) %>% summarise(n())
  plot_geu<-cbind(plot_geu,rr[,2])
}}
colnames(plot_geu)<-c("sample_id",">=8",">=16",">=30",">=60")
plot_geu_long<-plot_geu %>%
  pivot_longer(!sample_id, names_to = "group", values_to = "total")

for(i in c(8,16,30,60)){
  if(i==8){
    plot_SRA<-plot %>% group_by(indv) %>%  filter(total>=i) %>% summarise(n())}else{
      rr<-plot %>% group_by(indv) %>%  filter(total>=i) %>% summarise(n())
      plot_SRA<-cbind(plot_SRA,rr[,2])
    }}
colnames(plot_SRA)<-c("sample_id",">=8",">=16",">=30",">=60")
plot_SRA_long<-plot_SRA %>%
  pivot_longer(!sample_id, names_to = "group", values_to = "total")

plot_geu_long$study<-"Geuvadis"
plot_SRA_long$study<-"SRP187978"

join_plot<-rbind(plot_geu_long,plot_SRA_long)
level_order <- c('>=8', '>=16', '>=30',">=60")
pdf(file="~/plot/ASE/hetSNP.pdf", width = 10, height = 4)
ggplot()+
  geom_boxplot(data=join_plot, aes(x=group, y=total,color=study))+
  scale_x_discrete(limits = level_order)+
  labs(title="Comparison of total number of het_SNPs in 3 samples from Geuvadis vs SRP187978",
       subtitle="Runs for 3 individuals in SRA are combined",
       x="Reads/hetSNP",
       y="Number of het_SNP per indiv")
dev.off()
  

#----------------
#Get gene annotations
#------------------
library(AnnotationHub)
ah <- AnnotationHub()
ah <- query(ah, c("v26","GENCODE","Homo sapiens","GRch38")) #v26 was used in GTExV8
TxDb<-ah[["AH75155"]]
gene_grch38<-genes(TxDb,columns=c("TXID", "TXNAME"))


ase_gr<-makeGRangesFromDataFrame(geuvadis, seqnames.field="chr",
                                 start.field="pos", end.field="pos")

ov<-findOverlaps(gene_grch38,ase_gr)
geuvadis$gene_id<-NA
geuvadis$gene_id[subjectHits(ov)]<-names(gene_grch38)[queryHits(ov)]

geuvadis %>% group_by(sample_id) %>% summarise(length(unique(gene_id)))
geuvadis %>% group_by(gene_id, sample_id) %>%
  summarize(total=n()) %>% ungroup()


ase_gr<-makeGRangesFromDataFrame(plot, seqnames.field="chr",
                                 start.field="pos", end.field="pos")

ov<-findOverlaps(gene_grch38,ase_gr)
plot$gene_id<-NA
plot$gene_id[subjectHits(ov)]<-names(gene_grch38)[queryHits(ov)]

plot %>% summarise(length(unique(gene_id)))

for(i in c(8,16,30,60)){
  if(i==8){
    plot_geu<-geuvadis %>% group_by(sample_id) %>%  filter(total>=i) %>% summarise(length(unique(gene_id)))}else{
      rr<-geuvadis %>% group_by(sample_id) %>%  filter(total>=i) %>% summarise(length(unique(gene_id)))
      plot_geu<-cbind(plot_geu,rr[,2])
    }}
colnames(plot_geu)<-c("sample_id",">=8",">=16",">=30",">=60")
plot_geu_long<-plot_geu %>%
  pivot_longer(!sample_id, names_to = "group", values_to = "total")

for(i in c(8,16,30,60)){
  if(i==8){
    plot_SRA<-plot %>% group_by(indv) %>%  filter(total>=i) %>% summarise(length(unique(gene_id)))}else{
      rr<-plot %>% group_by(indv) %>%  filter(total>=i) %>% summarise(length(unique(gene_id)))
      plot_SRA<-cbind(plot_SRA,rr[,2])
    }}
colnames(plot_SRA)<-c("sample_id",">=8",">=16",">=30",">=60")
plot_SRA_long<-plot_SRA %>%
  pivot_longer(!sample_id, names_to = "group", values_to = "total")

plot_geu_long$study<-"Geuvadis"
plot_SRA_long$study<-"SRP187978"

join_plot<-rbind(plot_geu_long,plot_SRA_long)
level_order <- c('>=8', '>=16', '>=30',">=60")
pdf(file="~/plot/ASE/coding_genes.pdf", width = 10, height = 4)
ggplot()+
  geom_boxplot(data=join_plot, aes(x=group, y=total,color=study))+
  scale_x_discrete(limits = level_order) +
  labs(title="Comparison of total number of genes in 3 samples from Geuvadis vs SRP187978",
       subtitle="Runs for 3 individuals in SRA are combined",
       x="Reads/hetSNP",
       y="Number of genes")
dev.off()


#--------------
#cumsum plot
#--------------
colnames(plot)[3]<-"sample_id"
plot3<-plot %>% dplyr::select(c(sample_id,total))
plot3<-geuvadis %>% dplyr::select(c(sample_id,total)) %>% rbind(plot3)
pdf(file="~/plot/ASE/Total_cumulativeDis.pdf", width = 10, height = 4)
ggplot(plot3,aes(x=total,color=sample_id)) +
  stat_ecdf(geom="step")+
  xlim(c(0,200))+
  scale_color_manual(values=c("salmon","purple","pink","lightblue","steelblue","steelblue1"))+
  labs(title="Comparison of total read counts in 3 samples from Geuvadis vs SRP187978",
       subtitle="Runs for 3 individuals in SRA are combined",
       x="Reads/hetSNP",
       y="Cumulative distribution")
dev.off()

plot3<-plot %>%  rowwise() %>% mutate(min_allele=min(ref,alt)) %>% dplyr::select(c(sample_id,min_allele, total))
plot3<-geuvadis %>%rowwise() %>% mutate(min_allele=min(ref,alt)) %>% dplyr::select(c(sample_id,min_allele,total)) %>% rbind(plot3)

pdf(file="~/plot/ASE/minAllele_cumulativeDis.pdf", width = 10, height = 4)
ggplot(plot3,aes(x=min_allele,color=sample_id)) +
  stat_ecdf(geom="step")+
  xlim(c(0,200))+
  scale_color_manual(values=c("salmon","purple","pink","lightblue","steelblue","steelblue1"))+
  labs(title="Comparison of min_allele in 3 samples from Geuvadis vs SRP187978",
       subtitle="Runs for 3 individuals in SRA are combined",
       x="Reads/hetSNP",
       y="Cumulative distribution")
dev.off()

plot3.1<-plot3 %>% filter(total<=200)
cut(plot3$total,break=10)
pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
ggplot(plot3,aes(x=min_allele,y=total,color=sample_id)) +
  geom_line()+
  xlim(c(0,200))+
  ylim(c(0,1000))
dev.off()




