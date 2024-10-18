# In this analysis I am using the ASE data from GTEx
#This data was downloaded by Nick during rotations from dbGap: V8 not WASP corrected :location /dcl01/hansen/data/gtex_ase
#The GTEx paper recommends using the ase WASP corrected (allelic mapping error). Afrooz downloaded this on 08/09/2023 and is here /dcl01/hansen/data/arazi/ASE/dbGap

###For now start the analysis with ase not WASP corrected:

setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)
library(ggplot2)
library(cowplot)
library(scattermore)

#########################
#Plot data analysis: Liver
#########################
tissue<-"Liver"

ase_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", tissue, ".rds"))
true_gtex<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/true_gtex_", tissue, ".rds"))

colnames(true_gtex)[1:2]<-c("chr","pos")
#Inner join to only select the locations that are same between the two dataset
ase_join<-inner_join(true_gtex,ase_df, by=c("chr", "pos", "sample_id"))

#Take a look at how many snps and ase hits are in each data:
ase_join %>% group_by(sample_id) %>% summarize(ref_ratio=median(ref_ratio),n_snp=n(),
                                               recount_hit= sum(q_val<0.01,na.rm = T),
                                               gtex_hit=sum(BINOM_P_ADJUSTED<0.01,na.rm = T),
                                               same=sum(BINOM_P_ADJUSTED<0.01 & q_val<0.01,na.rm = T))
ase_sample<-ase_join %>% filter(sample_id==unique(ase_join$sample_id)[1])


x<-ase_join %>% filter(sample_id==unique(ase_df$sample_id)[1])%>% 
  rowwise() %>% mutate(mine_allele=min(alt,ref))

cov8<-data.frame(binom_p,binom_q)

# perform a binomial test for deviation from 0.5
binom_p = apply(x[,c("ref","alt")], 1, function(x) binom.test(x[1],x[1]+x[2],p=0.5)$p.value)
# perform multiple testing correction with FDR
binom_q = p.adjust(binom_p, method = "fdr")
x$cov8<-binom_q
sum(cov8$binom_q<0.01)
sum(x$q_val<0.01)
#---------------------------------------------
#Plot
#---------------------------------------------
pdf(file="~/plot/ASE/mono_allelic_filter.pdf", width = 10, height = 4)
for(i in unique(ase_df$sample_id)[1:3]){
  x<-ase_join %>% filter(sample_id==i)%>% rowwise() %>% mutate(mine_allele=min(alt,ref))

p=ggplot(x)+
  geom_point(aes(x=total, y=mine_allele))+xlim(0,5000)+ylim(0,3000)
print(p)
}
dev.off()
pdf(file="~/plot/ASE/test2.pdf", width = 10, height = 4)
ggplot(x)+
  geom_point(data=x %>% filter(q_val<0.01),aes(x=log10(ref), y=log10(alt)), color="red", alpha=0.4)+
  geom_point(data=x %>% filter(cov8<0.01),aes(x=log10(ref), y=log10(alt)), color="blue", alpha=0.4)
dev.off()
######
##Look at the p-value distribution across the samples:
#####
pdf(file="~/plot/ASE/p_q_val.pdf", width = 10, height = 4)
ggplot()+geom_histogram(data=ase_join, aes(x=p_val))+ facet_wrap(vars(sample_id), nrow = 3)+labs(title="Recount")
ggplot()+geom_histogram(data=ase_join, aes(x=q_val))+ facet_wrap(vars(sample_id), nrow = 3)

ggplot()+geom_histogram(data=ase_join, aes(x=BINOM_P))+ facet_wrap(vars(sample_id), nrow = 3)+labs(title="GTEx")
ggplot()+geom_histogram(data=ase_join, aes(x=BINOM_P_ADJUSTED))+ facet_wrap(vars(sample_id), nrow = 3)
dev.off()

######
##Look at the reference bias distribution across the samples:
#####
pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
ggplot()+
  geom_point(data=ase_join %>% filter(sample_id==sam),aes(x=-log2(BINOM_P_ADJUSTED), y= -log2(q_val)), alpha=0.5, color="red")+
  xlim(0,100)
dev.off()
pdf(file="~/plot/ASE/ref_ratio.pdf", width = 10, height = 4)
ggplot(ase_df) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+labs(title="all recount")

ggplot(ase_join) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+labs(title="recount")

ggplot(ase_join) + 
  geom_boxplot(aes(x=sample_id, y=REF_RATIO),position=position_dodge(1))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+labs(title="GTEx")

dev.off()



merge_gr_ase<-subset(merge_gr, BINOM_P_ADJUSTED<=0.05)
plot_df<-as.data.frame(merge_gr_ase)
pdf(file="~/plot/ASE/gtex_vs_recount.pdf", width = 6, height = 4)
ggplot()+
  geom_point(data=ase_join %>% filter(sample_id==sam,BINOM_P_ADJUSTED<0.01),aes(x=log10(REF_COUNT), y= log10(ALT_COUNT)), alpha=0.5, color="red")

ggplot()+
  geom_point(data=ase_join %>% filter(sample_id==sam,q_val<0.01), aes(x=log10(ref), y= log10(alt)), alpha=0.5, color="red")

ggplot()+
  geom_point(data=ase_join %>% filter(sample_id==unique(ase_df$sample_id)[2],BINOM_P_ADJUSTED<0.01),aes(x=log10(ref), y= log10(alt)), alpha=0.5, color="red")

ggplot()+
  geom_point(data=ase_join %>% filter(sample_id==unique(ase_df$sample_id)[2],q_val<0.01), aes(x=log10(REF_COUNT), y= log10(ALT_COUNT)), alpha=0.5, color="red")


dev.off()



#---------------------
#MA plot
#--------------------
t<-ase_join %>% 
  mutate(ratio=log2(TOTAL_COUNT) - log2(total),
          mean=(log2(TOTAL_COUNT) + log2(total))/2,
         ratio_ref=log2(REF_COUNT) - log2(ref),
         mean_ref=(log2(REF_COUNT) + log2(ref))/2,
         ratio_alt=log2(ALT_COUNT) - log2(alt),
         mean_alt=(log2(ALT_COUNT) + log2(alt))/2)

pdf(file="~/plot/ASE/MA_ref.pdf", width = 6, height = 4)
for (i in 1:10) {
  i
  p=ggplot(t %>% filter(sample_id == unique(ase_join$sample_id)[i]),aes(x=mean_ref, y= ratio_ref))+
    geom_point( alpha=0.1)+geom_smooth() +
    labs(x="log2(gtex+recount)/2", y="log2(gtex-recount)")
  print(p)
}
dev.off()

pdf(file="~/plot/ASE/MA_total.pdf", width = 6, height = 4)

  p=ggplot(t,aes(x=mean, y= ratio, color=sample_id),alpha=0.4)+
    geom_smooth() +
    labs(x="log2(gtex+recount)/2", y="log2(gtex-recount)")+ theme(legend.position='none')+
    labs(title = "Total MA plot of 47 first samples in Liver")
  print(p)

p=ggplot(t,aes(x=mean_ref, y= ratio_ref, color=sample_id),alpha=0.4)+
  geom_smooth() +
  labs(x="log2(gtex+recount)/2", y="log2(gtex-recount)")+ theme(legend.position='none')+
  labs(title = "REF MA plot of 47 first samples in Liver")
print(p)

p=ggplot(t,aes(x=mean_alt, y= ratio_alt, color=sample_id),alpha=0.4)+
  geom_smooth() +
  labs(x="log2(gtex+recount)/2", y="log2(gtex-recount)")+ theme(legend.position='none')+
  labs(title = "ALT MA plot of 47 first samples in Liver")
print(p)
dev.off()

#------------------------
#Loci analysis
#------------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene

m_gr<-makeGRangesFromDataFrame(ase_df, seqnames.field= "chr", start.field="pos", end.field="pos", ignore.strand=FALSE)
mm<-locateVariants(m_gr,txdb,AllVariants())
mm<-mm[,-c(2:10)]

m_gr<-makeGRangesFromDataFrame(ase_join_r, seqnames.field= "chr", start.field="pos", end.field="pos", ignore.strand=FALSE)
ov<-findOverlaps(m_gr, mm)

ase_join_r$location<-NA
ase_join_r$location[queryHits(ov)]<-as.character(mm$LOCATION[subjectHits(ov)])

pdf(file="~/plot/ASE/Gene_location.pdf", width = 10, height = 4)

ggplot(ase_join_r %>% filter(sample_id %in% unique(ase_join_r$sample_id)[1:5])) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio, color=location))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")

ggplot(ase_join_r %>% filter(sample_id %in% unique(ase_join_r$sample_id)[5:10])) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio, color=location))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")

ggplot(ase_join_r %>% filter(sample_id %in% unique(ase_join_r$sample_id)[10:15])) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio, color=location))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")

ggplot(ase_join_r %>% filter(sample_id %in% unique(ase_join_r$sample_id)[15:20])) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio, color=location))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")

ggplot(ase_join_r %>% filter(sample_id %in% unique(ase_join_r$sample_id)[20:25])) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio, color=location))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")

ggplot(ase_join_r %>% filter(sample_id %in% unique(ase_join_r$sample_id)[25:30])) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio, color=location))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")

ggplot(ase_join_r %>% filter(sample_id %in% unique(ase_join_r$sample_id)[30:35])) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio, color=location))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")

ggplot(ase_join_r %>% filter(sample_id %in% unique(ase_join_r$sample_id)[35:40])) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio, color=location))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")

ggplot(ase_join_r %>% filter(sample_id %in% unique(ase_join_r$sample_id)[40:47])) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio, color=location))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")
dev.off()

S<-ase_join_r %>% group_by(sample_id, location) %>% summarize(n=n())
S2<-ase_join_r%>% group_by(sample_id, VARIANT_ANNOTATION) %>% summarize(n=n())
S3<-ase_join %>% group_by(sample_id, VARIANT_ANNOTATION) %>% summarize(n=n())
pdf(file="~/plot/ASE/Gene_location_bar.pdf", width = 10, height = 4)


ggplot(data=S, aes(x=sample_id, y=n, fill=location)) +
  geom_bar(stat="identity")+ scale_x_discrete(labels=1:47)

ggplot(data=S2, aes(x=sample_id, y=n, fill=VARIANT_ANNOTATION)) +
  geom_bar(stat="identity")+ scale_x_discrete(labels=1:47)

ggplot(data=S3, aes(x=sample_id, y=n, fill=VARIANT_ANNOTATION)) +
  geom_bar(stat="identity")+ scale_x_discrete(labels=1:47)
dev.off()

sam<-unique(ase_join$sample_id)[29:33]

S2<-ase_join_r%>% filter(sample_id== sam)
S3<-ase_join %>%filter(sample_id== sam)
dim(S2)
dim(S3)
pdf(file="~/plot/ASE/abnormal.pdf", width = 10, height = 4)
ggplot(ase_join_r %>% filter(sample_id %in% sam)) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio, color=location))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  labs(title="Recount ref ratio")

ggplot(ase_join %>% filter(sample_id %in% sam)) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio, color=VARIANT_ANNOTATION))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  labs(title="gtex ref ratio")

ggplot(ase_join_r %>% filter(sample_id %in% sam)) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  labs(title="Recount ref ratio")

ggplot(ase_join %>% filter(sample_id %in% sam)) + 
  geom_boxplot(aes(x=sample_id, y=ref_ratio))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  labs(title="gtex ref ratio")

dev.off()
