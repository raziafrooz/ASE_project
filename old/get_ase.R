
setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)

#Gtex paper suggests removing vriants in HLA genes:
#his BED file contains genomic positions that we have identified as either showing bias in simulations or having a UCSC mappability score < 50. 
#Variants that fall into these positions are used for phasing, but not for generating haplotypic counts to avoid problems with mapping bias.
#https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
#Downloaded blacklist of snps from gtex (https://github.com/secastel/phaser/blob/master/phaser/README.md)
bad_snp<-read.table("~/plot/ASE/hg38_haplo_count_blacklist.chr.bed", sep="\t")
bad_snp_gr<-makeGRangesFromDataFrame(bad_snp,seqnames="V1",start.field ="V2",end.field = "V3")


tissues<-c("Liver", "Lung", "Stomach", "Pancreas")

for(k in 1:length(tissues)){
  tissue<-tissues[k]
  if(!file.exists(file=paste0("~/hansen_lab/ASE/test_ASE/", tissue, "_ase.rda") ))
  {
  print(tissue)
ase_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", tissue, ".rds"))
ase_df$ref_ratio<- ase_df$ref/ase_df$total
for(i in 1:length(unique(ase_df$sample_id))){
  print(i)
  sam<-unique(ase_df$sample_id)[i]
  
  ase_filt<-ase_df %>% filter(sample_id==sam)
  ase_filt_gr<-makeGRangesFromDataFrame(ase_filt,seqnames="chr",start.field ="pos",end.field = "pos", keep.extra.columns = T)
  
#Remove blacklist snps from our granges:
ov<-findOverlaps(ase_filt_gr,bad_snp_gr)
ase_filt_gr<-ase_filt_gr[-unique(queryHits(ov))]

#obtain the median of the ref ratio to be used as the null p value
is.median<-median(ase_filt_gr$ref_ratio)
ase_filt_gr$p_val = apply(as.data.frame(ase_filt_gr)[,c("ref","alt")], 1, function(x) binom.test(x[1],(x[1]+x[2]),p=is.median)$p.value)
# perform multiple testing correction with FDR
ase_filt_gr$q_val = p.adjust(ase_filt_gr$p_val, method = "fdr")

if(i==1){
  ase_all<- as.data.frame(ase_filt_gr)
} else{
  ase_all<-rbind(ase_all,as.data.frame(ase_filt_gr))
}
}
#order the df and save the data:
colnames(ase_all)[1:2]<-c("chr","pos")
ase_all[,3:5]<-NULL
ase_all<-as_tibble(ase_all)
save(ase_all, file=paste0("~/hansen_lab/ASE/test_ASE/", tissue, "_ase.rda") )
  }}

k=3
#check ref ratio plots
pdf(file="~/plot/ASE/ref_ratio.pdf", width = 10, height = 4)
for(k in 1:4){
  tissue<-tissues[k]
  print(tissue)
load(paste0("~/hansen_lab/ASE/test_ASE/", tissue, "_ase.rda"))
true_gtex<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/true_gtex_", tissue, ".rds"))
colnames(true_gtex)[1:2]<-c("chr","pos")

#ase_join_filter<-ase_join
#compare the results in recount and gtex:
ase_join<-inner_join(true_gtex,ase_all, by=c("chr", "pos", "sample_id"))
#x2<-ase_join %>% group_by(sample_id) %>% summarize(med=median(ref_ratio),
#                                                  med_gtex=median(REF_RATIO))
#x2$tissue<-tissue
#x<-rbind(x,x2)}
p<-ggplot(ase_join)+
  geom_boxplot(aes(x=sample_id,y=REF_RATIO),alpha=0.02)+ 
  labs(title=paste0("ref ratio in GTEx: ", tissue))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")
print(p)

p<-ggplot(ase_join)+
  geom_boxplot(aes(x=sample_id,y=ref_ratio),alpha=0.02)+ 
  labs(title=paste0("ref ratio in Recount: ", tissue))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")
print(p)
}
dev.off()
ase_join
x<-ase_join %>% group_by(sample_id) %>% summarize(med=median(ref_ratio),
                                               med_gtex=median(REF_RATIO))
pdf(file="~/plot/ASE/ref_ratio_all.pdf", width = 10, height = 4)
ggplot(x,aes(x=med_gtex,y=med, color=tissue))+
  geom_point(alpha=0.5)+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  geom_vline(xintercept = 0.5, color="red",linetype="dotdash")
dev.off()
summary(x$med_gtex)
x <- ase_join %>% filter(sample_id %in% unique(ase_all$sample_id))
x2<-x%>% filter(q_val<0.01)%>% select("chr", "pos", "sample_id", "q_val") %>% pivot_wider(names_from= sample_id,
                  values_from = q_val)

x3<-x2 %>% summarize(rowSums(is.na(.)))
x3<-x2 %>% drop_na()
str(x3)
id<-which(x3<27)
x2[id,]

#---------------------------------------
#Get gene annotations
#---------------------------------------
#library(AnnotationHub)
#ah <- AnnotationHub(localHub=FALSE)

# query(ah, "Gencode", "v29")
# gencode_29<-ah[["AH75174"]]
x <- ase_join %>% filter(sample_id %in% unique(ase_all$sample_id)[1:3])
x2<-x %>% group_by(sample_id, GENE_ID) %>% summarize(n=n(),
                                                 hit_total=sum(q_val<0.01),
                                                 hit_total_gtex=sum(BINOM_P_ADJUSTED<0.01))

pdf(file="~/plot/ASE/gene.pdf", width = 10, height = 4)
ggplot(x2)+
  geom_bar(aes(n))+facet_wrap(vars(sample_id))
dev.off()
tail(x2)
x2 %>% filter(hit_total>0)
quantile()
ase_join$total_cut[which(ase_join$total>=8)]<-">=8"
ase_join$total_cut[which(ase_join$total>=16)]<-">=16"
ase_join$total_cut[which(ase_join$total>=30)]<-">=30"
ase_join$total_cut[which(ase_join$total>=60)]<-">=60"
ase_join$total_cut<-factor(ase_join$total_cut, levels = c(">=8",">=16",">=30",">=60"))
ase_join<-ase_join %>% group_by(sample_id,total_cut) %>%  mutate(n_het_snp=n())%>% ungroup()

ase_join$total_cut_g[which(ase_join$TOTAL_COUNT>=8)]<-">=8"
ase_join$total_cut_g[which(ase_join$TOTAL_COUNT>=16)]<-">=16"
ase_join$total_cut_g[which(ase_join$TOTAL_COUNT>=30)]<-">=30"
ase_join$total_cut_g[which(ase_join$TOTAL_COUNT>=60)]<-">=60"
ase_join$total_cut_g<-factor(ase_join$total_cut_g, levels = c(">=8",">=16",">=30",">=60"))
ase_join<-ase_join %>% group_by(sample_id,total_cut_g) %>%  mutate(n_het_snp_gtex=n()) %>% ungroup()

pdf(file="~/plot/ASE/cov3.1.pdf", width = 10, height = 4)
ggplot(ase_join)+
  geom_boxplot(aes(x=total_cut, y=n_het_snp))
ggplot(ase_join)+
  geom_boxplot(aes(x=total_cut_g, y=n_het_snp_gtex))

dev.off()

pdf(file="~/plot/ASE/mono_allelic_filter.pdf", width = 10, height = 4)
for(i in unique(ase_df$sample_id)[1:3]){
  x<-ase_join %>% filter(sample_id==i)%>% rowwise() %>% mutate(mine_allele=min(alt,ref))
  
  p=ggplot(x)+
    geom_point(aes(x=total, y=mine_allele))+xlim(0,5000)+ylim(0,3000)
  print(p)
}
dev.off()



####---------------------
#QQplot P-val comparison
####---------------------
# ggd.qqplot = function(pvector, main=NULL, ...) {
#   o = -log10(sort(pvector,decreasing=F))
#   e = -log10( 1:length(o)/length(o) )
#   plot(e,o,pch=19,cex=1, main=main, ...,
#        xlab=expression(Expected~~-log[10](italic(p))),
#        ylab=expression(Observed~~-log[10](italic(p))),
#        xlim=c(0,6), ylim=c(0,60))
#   lines(e,e,col="red")
# }
o = -log10(sort(ase_join_filter$p_val,decreasing=F))
e = -log10( 1:length(o)/length(o) )
o2=-log10(sort(ase_join$p_val,decreasing=F))
e2=-log10( 1:length(sort(ase_join$p_val,decreasing=F))/length(sort(ase_join$p_val,decreasing=F)))
pdf(file="~/plot/ASE/qqplot.pdf", width = 10, height = 4)
plot(e,o,pch=19,cex=0.4, main="p-val QQ plot",
     xlab=expression(Expected~~-log[10](italic(p))),
     ylab=expression(Observed~~-log[10](italic(p))),
     xlim=c(0,4), ylim=c(0,60),col=alpha("darkgreen",0.5))
points(e2,o2,pch=19,cex=0.4,col=alpha("lightblue",0.5))
lines(e,e,col="red")
dev.off()


#---------------------------
#calculate the AE and plot GTEx vs recount3
#---------------------------
x<-ase_join %>% filter(sample_id == unique(ase_join$sample_id)[1])
ase_join$ae<- abs(0.5-(ase_join$ref/ase_join$total))
ase_join$AE<- abs(0.5-ase_join$REF_RATIO)
pdf(file="~/plot/ASE/AE_comparison.pdf", width = 10, height = 4)
ggplot(ase_join)+
  geom_point(aes(x= AE, y=ae))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")

ggplot(ase_join %>% filter(BINOM_P_ADJUSTED<0.01))+
  geom_point(aes(x= AE, y=ae))+
   geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")

ggplot(ase_join %>% filter(p_val<0.01))+
  geom_point(aes(x= AE, y=ae))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")

dev.off()

pdf(file="~/plot/ASE/ref_ratio_comparison.pdf", width = 10, height = 4)
ggplot(ase_join)+
  geom_point(aes(x= REF_RATIO, y=(ref/total)),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")

ggplot(ase_join %>% filter(BINOM_P_ADJUSTED<0.01))+
  geom_point(aes(x= REF_RATIO, y=(ref/total)),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")

ggplot(ase_join %>% filter(p_val<0.01))+
  geom_point(aes(x= REF_RATIO, y=(ref/total)),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")

dev.off()
n<-as.numeric(quantile(ase_join$total))
ase_join$total_cut<-cut(ase_join$total, breaks=c(n[1],n[2],n[3],n[4], n[5]), include.lowest=T)

pdf(file="~/plot/ASE/AE_total_cut.pdf", width = 10, height = 4)
ggplot(ase_join)+
  geom_point(aes(x= AE, y=ae),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+ 
  facet_wrap(vars(total_cut))

ggplot(ase_join %>% filter(BINOM_P_ADJUSTED<0.01))+
  geom_point(aes(x= AE, y=ae),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))

ggplot(ase_join %>% filter(p_val<0.01))+
  geom_point(aes(x= AE, y=ae),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))

dev.off()
pdf(file="~/plot/ASE/refRatio_total_cut.pdf", width = 10, height = 4)
ggplot(ase_join)+
  geom_point(aes(x= REF_RATIO, y=(ref/total)),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+ 
  facet_wrap(vars(total_cut))

ggplot(ase_join %>% filter(BINOM_P_ADJUSTED<0.01))+
  geom_point(aes(x= REF_RATIO, y=(ref/total)),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))

ggplot(ase_join %>% filter(p_val<0.01))+
  geom_point(aes(x= REF_RATIO, y=(ref/total)),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))

dev.off()

n<-as.numeric(quantile(ase_join$TOTAL_COUNT))
ase_join$TOTAL_cut<-cut(ase_join$TOTAL_COUNT, breaks=c(n[1],n[2],n[3],n[4], n[5]), include.lowest=T)

pdf(file="~/plot/ASE/AE_TOTAL_cut.pdf", width = 10, height = 4)
ggplot(ase_join)+
  geom_point(aes(x= AE, y=ae),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+ 
  facet_wrap(vars(TOTAL_cut))

ggplot(ase_join %>% filter(BINOM_P_ADJUSTED<0.01))+
  geom_point(aes(x= AE, y=ae),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(TOTAL_cut))

ggplot(ase_join %>% filter(p_val<0.01))+
  geom_point(aes(x= AE, y=ae),alpha=0.4)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(TOTAL_cut))

dev.off()






         