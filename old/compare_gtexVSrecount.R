setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)
tissue<-"Liver"

load(paste0("~/hansen_lab/ASE/test_ASE/", tissue, "_ase.rda") )
true_gtex<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/true_gtex_", tissue, ".rds"))
colnames(true_gtex)[1:2]<-c("chr","pos")

sam<- unique(ase_all$sample_id)[2]
ase_all_1 <- ase_all %>% filter(sample_id %in% sam)
true_gtex_1 <- true_gtex %>% filter(sample_id %in% sam)

# calculate AE:
ase_all_1$ae<-abs(0.5-ase_all_1$ref_ratio)
true_gtex_1$AE<-abs(0.5-true_gtex_1$REF_RATIO)

ase_all_1$ae<-abs((ase_all_1$alt + 1)/(ase_all_1$ref + 1))
true_gtex_1$AE<-abs((true_gtex_1$ALT_COUNT + 1)/(true_gtex_1$REF_COUNT + 1))

ase_all_1$aFC_recount<-log2((ase_all_1$alt + 1)/(ase_all_1$ref + 1))
true_gtex_1$aFC_gtex<-log2((true_gtex_1$ALT_COUNT + 1)/(true_gtex_1$REF_COUNT + 1))


#pick the SNPs that are the same:
ase_same<-inner_join(true_gtex_1,ase_all_1, by=c("chr", "pos", "sample_id"))
ase_same$ae_difference<-ase_same$AE-ase_same$ae
#separate based on the coverage
n<-as.numeric(quantile(ase_same$TOTAL_COUNT))
ase_same$TOTAL_cut<-cut(ase_same$TOTAL_COUNT, breaks=c(n[1],n[2],n[3],n[4], n[5]), include.lowest=T)


#--------------------------------------
#Make plots:
#--------------------------------------
x<-ase_same %>% group_by(sample_id) %>% rowwise() %>% mutate(recount_min=min(ref,alt), gtex_min= min(REF_COUNT,ALT_COUNT),
                                     ASE_sig=case_when(BINOM_P<0.05 ~ "sig", BINOM_P>0.05 ~ "insig"),
                                     ratio=(aFC_gtex) - (aFC_recount),
                                     mean=((aFC_gtex) + (aFC_recount))/2,
                                     ratio2=log2(AE) - log2(ae),
                                     mean2=(log2(AE) + log2(ae))/2)

n<-as.numeric(quantile(x$ratio))
x$aFC_cut<-cut(x$ratio, breaks=c(n[1],n[2],n[3],n[4], n[5]), include.lowest=T)
x$dif<-abs(x$recount_min-x$gtex_min)

x$ref_dif<-x$REF_COUNT-x$ref
x$alt_dif<-x$ALT_COUNT-x$alt
x$ae_dif<-x$AE-x$ae
pdf(file="~/plot/ASE/aFC_ratio_cut_multiSample.pdf", width = 10, height = 6)
p=ggplot(x,aes(x=TOTAL_COUNT,y=ratio, color=sample_id))+
  theme(legend.position="none")+
  stat_smooth()+
  lims(x=c(0,50), y=c(-0.75,0.75))+
  labs(title="All SNPs multiple samples", x="gtex total count", y="aFC Log ratio")
print(p)

p=ggplot(x%>% filter(BINOM_P<0.05),aes(x=alt_dif,y=ref_dif, color=sample_id))+
  theme(legend.position="none")+
  stat_smooth()+
  #geom_point(alpha=0.2)+
  facet_wrap(vars(aFC_cut))+
  lims(x=c(-25,25), y=c(-50,50))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  labs(title="Gtex sig SNPs. 30 samples", subtitle="facet based on the logratio quantiles", x="gtex(alt) - recount(alt)", y="gtex(ref) - recount(ref)")
print(p)

p=ggplot(x%>% filter(BINOM_P>0.05),aes(x=alt_dif,y=ref_dif, color=sample_id))+
  theme(legend.position="none")+
  stat_smooth()+
  #geom_point(alpha=0.2)+
  facet_wrap(vars(aFC_cut))+
  lims(x=c(-25,25), y=c(-50,50))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  labs(title="Gtex insig SNPs. 30 samples", subtitle="facet based on the logratio quantiles", x="gtex(alt) - recount(alt)", y="gtex(ref) - recount(ref)")
print(p)
dev.off()
pdf(file="~/plot/ASE/aFC_mean_ratio.pdf", width = 10, height = 6)
ggplot(x)+
  geom_point(aes(x=TOTAL_COUNT, y=ratio),alpha=0.4)+
  facet_wrap(vars(ASE_sig))+lims(x=c(0,2000),y=c(-2,2))+
  labs(title="One sample all SNPs: aFC vs total count in Gtex")
dev.off()
pdf(file="~/plot/ASE/aFC_countDifference.pdf", width = 10, height = 6)
ggplot(x %>%  filter(sample_id== "GTEX-S33H-1626-SM-4AD68.1"))+
  geom_point(aes(x=TOTAL_COUNT, y=ratio),alpha=0.4)+
  facet_wrap(vars(ASE_sig))+lims(x=c(0,2000),y=c(-2,2))+
  labs(title="One sample all SNPs: aFC vs total count in Gtex")

ggplot(x %>% filter(ASE_sig=="sig",sample_id== "GTEX-S33H-1626-SM-4AD68.1"),aes(x=dif, y=ratio))+
  geom_point(alpha=0.4)+
  facet_wrap(vars(TOTAL_cut))+lims(x=c(0,20),y=c(-2,2))+
  geom_quantile(color="red")+
  labs(title="One sample sig SNPs: aFC vs difference of alelle with lower count between Recount and Gtex",
       subtitle="Facet based on gtex total quantiles\nRed lines are 20%-mean_75%")

ggplot(x %>% filter(ASE_sig=="insig",sample_id== "GTEX-S33H-1626-SM-4AD68.1"),aes(x=dif, y=ratio))+
  geom_point(alpha=0.4)+
  facet_wrap(vars(TOTAL_cut))+lims(x=c(0,20),y=c(-2,2))+
  geom_quantile( color="red")+
  labs(title="One sample insig SNPs: aFC vs difference of alelle with lower count between Recount and Gtex",
       subtitle="Facet based on gtex total quantiles\nRed lines are 20%-mean_75%")

ggplot(x %>% filter(ASE_sig=="sig"),aes(x=dif, y=ratio))+
  facet_wrap(vars(TOTAL_cut))+lims(x=c(0,20),y=c(-2,2))+
  geom_quantile(aes(color=sample_id))+
  theme(legend.position = "none")+
  labs(title="20 sample in Liver, sig SNPs: aFC vs difference of alelle with lower count between Recount and Gtex",
       subtitle="Facet based on gtex total quantiles\nRed lines are 20%-mean_75%")

ggplot(x %>% filter(ASE_sig=="insig"),aes(x=dif, y=ratio))+
  facet_wrap(vars(TOTAL_cut))+lims(x=c(0,20),y=c(-2,2))+
  geom_quantile(aes(color=sample_id))+
  theme(legend.position = "none")+
  labs(title="20 sample in Liver, insig SNPs: aFC vs difference of alelle with lower count between Recount and Gtex",
       subtitle="Facet based on gtex total quantiles\nRed lines are 20%-mean_75%")

dev.off()
ase_same<-ase_same %>% group_by(TOTAL_cut, sample_id) %>%  
  mutate(ratio=(aFC_gtex) - (aFC_recount),
         mean=((aFC_gtex) + (aFC_recount))/2,
         sd=sd(ratio))
x<-ase_same %>% mutate(ASE_sig=case_when(BINOM_P<0.05 ~ "sig",
         BINOM_P>0.05 ~ "insig")) %>% 
  group_by(TOTAL_cut,sample_id,ASE_sig) %>% summarize(v=var(ratio), median_ratio=median(ratio))

pdf(file="~/plot/ASE/aFC_mut.pdf", width = 10, height = 6)
for(i in 1: length(unique(ase_same$sample_id))){
  sam=unique(ase_same$sample_id)[i]
  p=ggplot(ase_same %>% filter(sample_id==sam,BINOM_P<0.05),aes(x=mean,y=ratio))+
    geom_point(alpha=0.4)+
  facet_wrap(vars(VARIANT_ANNOTATION))+
  geom_smooth()+
    labs(title=paste0("aFC MA plot: sig SNPs in ", sam))
print(p)

p=ggplot(ase_same %>% filter(sample_id==sam,BINOM_P>0.05),aes(x=mean,y=ratio))+
  geom_point(alpha=0.4)+
  facet_wrap(vars(VARIANT_ANNOTATION))+
  geom_smooth()+
  labs(title=paste0("aFC MA plot: insig SNPs in ", sam))
print(p)
}
dev.off()

x$ref_dif<-x$REF_COUNT-x$ref
x$alt_dif<-x$ALT_COUNT-x$alt

pdf(file="~/plot/ASE/aFC_ratio_cut2.pdf", width = 10, height = 6)

p=ggplot(x %>% filter(REF_COUNT>ALT_COUNT),aes(x=alt_dif,y=ref_dif))+
  geom_point(alpha=0.2)+
  geom_point(data=x %>% filter(REF_COUNT<ALT_COUNT),aes(x=alt_dif,y=ref_dif), color="blue",alpha=0.2)+
  facet_wrap(vars(aFC_cut))+
  lims(x=c(-50,50), y=c(-50,50))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  labs(title="All SNPs", subtitle="facet based on the logratio quantiles", x="gtex(alt) - recount(alt)", y="gtex(ref) - recount(ref)")
print(p)

p=ggplot(x%>% filter(BINOM_P<0.05,REF_COUNT>ALT_COUNT),aes(x=alt_dif,y=ref_dif))+
  geom_point(alpha=0.2)+
  geom_point(data=x %>% filter(BINOM_P<0.05, REF_COUNT<ALT_COUNT),aes(x=alt_dif,y=ref_dif), color="blue",alpha=0.2)+
  facet_wrap(vars(aFC_cut))+
  lims(x=c(-50,50), y=c(-50,50))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  labs(title="GTEx significant SNPs", subtitle="facet based on the logratio quantiles", x="gtex(alt) - recount(alt)", y="gtex(ref) - recount(ref)")
print(p)

p=ggplot(x%>% filter(BINOM_P>0.05,REF_COUNT>ALT_COUNT),aes(x=alt_dif,y=ref_dif))+
  geom_point(alpha=0.2)+
  geom_point(data=x %>% filter(BINOM_P>0.05, REF_COUNT<ALT_COUNT),aes(x=alt_dif,y=ref_dif), color="blue",alpha=0.2)+
  facet_wrap(vars(aFC_cut))+
  lims(x=c(-50,50), y=c(-50,50))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  labs(title="GTEx insignificant SNPs", subtitle="facet based on the logratio quantiles", x="gtex(alt) - recount(alt)", y="gtex(ref) - recount(ref)")
print(p)
dev.off()
pdf(file="~/plot/ASE/aFC_differences.2.pdf", width = 10, height = 6)

p=ggplot(x,aes(x=ae_dif,y=ratio2))+
  geom_point(alpha=0.2)+
  facet_wrap(vars(ASE_sig))+
  lims(x=c(-1,1))+
  labs(title="All SNPs", subtitle = "facet based on significant vs insignificant", x="gtex(ref) - recount(ref)", y="log2(gtex(alt+1/ref+1) - log2(recount(alt+1/ref+1))")
print(p)

p=ggplot(x %>% filter(BINOM_P<0.05, REF_COUNT>ALT_COUNT),aes(x=ae_dif,y=ratio2))+
  geom_point(alpha=0.2)+
  geom_point(data=x %>% filter(BINOM_P<0.05, REF_COUNT<ALT_COUNT),aes(x=ae_dif,y=ratio2), color="blue",alpha=0.2)+
  facet_wrap(vars(TOTAL_cut))+
  lims(x=c(-1,1))+
  labs(title="GTEx significant SNPs",subtitle="facet based on gtex total count. Blue dots are where gtex(ref<alt) count", x="gtex(ref) - recount(ref)", y="log2(gtex(alt+1/ref+1) - log2(recount(alt+1/ref+1))")
print(p)

p=ggplot(x %>% filter(BINOM_P>0.05,REF_COUNT>ALT_COUNT),aes(x=ae_dif,y=ratio2))+
  geom_point(alpha=0.2)+
  geom_point(data=x %>% filter(BINOM_P>0.05, REF_COUNT<ALT_COUNT),aes(x=ae_dif,y=ratio2), color="blue",alpha=0.2)+
  facet_wrap(vars(TOTAL_cut))+
  lims(x=c(-1,1))+
  labs(title="GTEx insignificant SNPs",subtitle="facet based on gtex total count. Blue dots are where gtex(ref<alt) count", x="gtex(ref) - recount(ref)", y="log2(gtex(alt+1/ref+1) - log2(recount(alt+1/ref+1))")
print(p)
dev.off()

p=ggplot(ase_same %>% filter(sample_id==sam,BINOM_P<0.05),aes(x=log2(ref),y=log2(alt)))+
  geom_point(alpha=0.2)+
  geom_point(aes(x=log2(REF_COUNT),y=log2(ALT_COUNT)), color="red",alpha=0.2)
print(p)
p=ggplot(ase_same %>% filter(sample_id==sam,BINOM_P<0.05),aes(x=log2(ref),y=log2(REF_COUNT)))+
  geom_point(alpha=0.2)
print(p)
dev.off()
pdf(file="~/plot/ASE/aFC.pdf", width = 10, height = 6)
for(i in 1: length(unique(ase_same$sample_id))){
  sam=unique(ase_same$sample_id)[i]
p=ggplot(ase_same %>% filter(sample_id==sam),aes(x=mean,y=ratio))+
  geom_point()+
  facet_wrap(vars(TOTAL_cut))+
    geom_smooth()+
    labs(title=paste0("aFC MA plot: all sample=", sam))
print(p)  
p= ggplot(ase_same %>% filter(sample_id==sam, BINOM_P<0.05),aes(x=mean,y=ratio))+
    geom_point()+
    facet_wrap(vars(TOTAL_cut))+
    geom_smooth()+
    labs(title=paste0("aFC MA plot: sig ASE sample=", sam))
  print(p)
  p=ggplot(ase_same %>% filter(sample_id==sam, BINOM_P>0.05),aes(x=mean,y=ratio))+
    geom_point()+
    facet_wrap(vars(TOTAL_cut))+
    geom_smooth()+
    labs(title=paste0("aFC MA plot: insig ASE sample=", sam))
  print(p)
  }
dev.off()


pdf(file="~/plot/ASE/aFC2.pdf", width = 10, height = 6)

ggplot(x %>% filter(ASE_sig =="sig", sample_id !="GTEX-ZF29-2026-SM-4WWB7.1"),aes(x=sample_id,y=median_ratio,ymin = median_ratio - v, ymax = median_ratio + v, color=TOTAL_cut))+
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(position = position_dodge(width = 0.5))+
  labs(title="aFC mean and variations plot: GTEx sig ASE")

ggplot(x%>% filter(ASE_sig =="insig",sample_id !="GTEX-ZF29-2026-SM-4WWB7.1"),aes(x=sample_id,y=median_ratio,ymin = median_ratio - v, ymax = median_ratio + v, color=TOTAL_cut))+
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(position = position_dodge(width = 0.5))+
  labs(title="aFC mean and variations plot: GTEx sig ASE")
dev.off()


#--------------------
n<-as.numeric(quantile(ase_same$aFC_recount))
ase_same$aFC_cut<-cut(ase_same$aFC_gtex, breaks=c(n[1],n[2],n[3],n[4], n[5]), include.lowest=T)


x<-ase_same %>% group_by(aFC_cut, sample_id) %>%  
  mutate(ratio=log2((TOTAL_COUNT)) - log2((total)),
         mean=(log2(TOTAL_COUNT) + log2(total))/2,
         ratio_ref=log2(REF_COUNT) - log2(ref),
         mean_ref=(log2(REF_COUNT) + log2(ref))/2,
         ASE_sig=case_when(BINOM_P<0.05 ~ "sig",
                                  BINOM_P>0.05 ~ "insig"))


pdf(file="~/plot/ASE/aFC_cut.pdf", width = 10, height = 6)

ggplot(x %>% filter(ASE_sig== "sig"),
       aes(x=mean, y=ratio, color=sample_id))+
  facet_wrap(vars(aFC_cut))+
  geom_smooth()+
  theme(legend.position = "none")
ggplot(x %>% filter(ASE_sig== "insig"),
       aes(x=mean, y=ratio, color=sample_id))+
  facet_wrap(vars(aFC_cut))+
  geom_smooth()+
  theme(legend.position = "none")

ggplot(x %>% filter(ASE_sig== "sig"),
       aes(x=mean_ref, y=ratio_ref, color=sample_id))+
  facet_wrap(vars(aFC_cut))+
  geom_smooth()+
  theme(legend.position = "none")

ggplot(x %>% filter(ASE_sig== "insig"),
       aes(x=mean_ref, y=ratio_ref, color=sample_id))+
  facet_wrap(vars(aFC_cut))+
  geom_smooth()+
  theme(legend.position = "none")
dev.off()

x<-ase_same %>% group_by(sample_id) %>%  
  mutate(ratio=log2((aFC_gtex)) - log2((aFC_recount)),
         mean=(log2(aFC_gtex) + log2(aFC_recount))/2,
         ASE_sig=case_when(BINOM_P<0.05 ~ "sig",
                           BINOM_P>0.05 ~ "insig"))


pdf(file="~/plot/ASE/aFC_MA.pdf", width = 10, height = 6)

ggplot(x %>% filter(sample_id==sam, ASE_sig== "sig"),
       aes(x=mean, y=ratio))+
  geom_point()+
  geom_smooth()
ggplot(x %>% filter(sample_id==sam, ASE_sig== "insig"),
       aes(x=mean, y=ratio))+
  geom_point()+
  geom_smooth()
dev.off()


ase_same %>% group_by(sample_id) %>% summarize(n=sum(BINOM_P<0.05),
                                               n_re=sum(p_val<0.05),
                                               n_all=sum(BINOM_P<0.05 & p_val<0.05))








ase_same<-ase_same %>% 
  mutate(ratio=(AE) - (ae),
         mean=((AE) + (ae))/2)

pdf(file="~/plot/ASE/allele_ratio.pdf", width = 10, height = 6)
ggplot(ase_same %>% filter(BINOM_P<0.05),aes(x= mean, y=ratio))+
  geom_point()+
  geom_smooth()+
  facet_wrap(vars(total_cut))+
  labs(title="Allele ratio in one sample: gtex significant",
       x="Mean", y="ratio (gtex_AE/recount_ae)")

ggplot(ase_same %>% filter(BINOM_P>0.05),aes(x= mean, y=ratio))+
  geom_point()+
  geom_smooth()+
  facet_wrap(vars(total_cut))+
  labs(title="Allele ratio in one sample: gtex insignificant",
       x="Mean", y="ratio (gtex_AE/recount_ae)")
dev.off()

pdf(file="~/plot/ASE/allele_ratio_recount.pdf", width = 10, height = 6)
ggplot(ase_same %>% filter(p_val<0.05),aes(x= mean, y=ratio))+
  geom_point()+
  geom_smooth()+
  facet_wrap(vars(total_cut))+
  labs(title="Allele ratio in one sample: recount significant",
       x="AE(GTEx)", y="ae(recount)")

ggplot(ase_same %>% filter(p_val>0.05),aes(x= mean, y=ratio))+
  geom_point()+
  geom_smooth()+
  facet_wrap(vars(total_cut))+
  labs(title="Allele ratio in one sample: recount insignificant",
       x="AE(GTEx)", y="ae(recount)")
dev.off()

plot_df<-ase_same %>% select(AE,ae,total_cut,BINOM_P,p_val) %>% pivot_longer(!c(BINOM_P,p_val,total_cut),names_to="dataset", values_to = "ratio")
plot_df$dataset[plot_df$dataset=="ae"]<-"recount"
plot_df$dataset[plot_df$dataset=="AE"]<-"GTEx"
pdf(file="~/plot/ASE/allele_ratio_box.pdf", width = 10, height = 6)
ggplot(plot_df %>% filter(BINOM_P<0.05))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))+
  labs(title="Allele ratio in one sample: GTEx significant",
       x="recount_total_cut", y="allelic ratio")

ggplot(plot_df %>% filter(BINOM_P>0.05))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))+
  labs(title="Allele ratio in one sample: GTEx insignificant",
       x="AE(GTEx)", y="ae(recount)")

ggplot(plot_df %>% filter(p_val<0.05))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))+
  labs(title="Allele ratio in one sample: recount significant",
       x="recount_total_cut", y="allelic ratio")

ggplot(plot_df %>% filter(p_val>0.05))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))+
  labs(title="Allele ratio in one sample: recount insignificant",
       x="recount_total_cut", y="allelic ratio")

ggplot(plot_df %>% filter(BINOM_P<0.05,p_val<0.05))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))+
  labs(title="Allele ratio in one sample: both significant",
       x="recount_total_cut", y="allelic ratio")

ggplot(plot_df %>% filter(BINOM_P>0.05,p_val>0.05))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))+
  labs(title="Allele ratio in one sample: both insignificant",
       x="recount_total_cut", y="allelic ratio")

ggplot(plot_df %>% filter(BINOM_P<0.05,p_val>0.05))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))+
  labs(title="Allele ratio in one sample: GTEx sig, recount insig",
       x="recount_total_cut", y="allelic ratio")

ggplot(plot_df %>% filter(BINOM_P>0.05,p_val<0.05))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))+
  labs(title="Allele ratio in one sample: GTEx insig, recount sig",
       x="recount_total_cut", y="allelic ratio")

dev.off()


pdf(file="~/plot/ASE/allele_ratio_dif.pdf", width = 10, height = 6)
ggplot(ase_same %>% filter(BINOM_P<0.05,p_val<0.05))+
  geom_point(aes(x= AE, y=ae))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))+
  labs(title="Allele ratio in one sample: both significant",
       x="AE(GTEx)", y="ae(recount)")

ggplot(ase_same %>% filter(BINOM_P<0.05,p_val>0.05))+
  geom_point(aes(x= AE, y=ae))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))+
  labs(title="Allele ratio in one sample: GTEx sig, recount insig",
       x="AE(GTEx)", y="ae(recount)")

ggplot(ase_same %>% filter(BINOM_P>0.05,p_val<0.05))+
  geom_point(aes(x= AE, y=ae))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))+
  labs(title="Allele ratio in one sample: GTEx insig, recount sig",
       x="AE(GTEx)", y="ae(recount)")
dev.off()


pdf(file="~/plot/ASE/ref_ratio_compare.pdf", width = 10, height = 6)
ggplot(ase_same)+
  geom_point(aes(x= REF_RATIO, y=ref_ratio,color=total_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))+
  labs(title="similar snps:all snps separated by recount total")

ggplot(ase_same %>% filter(BINOM_P<0.05 & p_val<0.05))+
  geom_point(aes(x= REF_RATIO, y=ref_ratio,color=total_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))+
  labs(title="similar snps:hit snps in both dataset (p<0.05) separated by recount total")

ggplot(ase_same %>% filter(BINOM_P<0.05))+
  geom_point(aes(x= REF_RATIO, y=ref_ratio,color=total_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))+
  labs(title="similar snps:gtex hit snps (p<0.05) separated by recount total")

ggplot(ase_same %>% filter(p_val<0.05))+
  geom_point(aes(x= REF_RATIO, y=ref_ratio,color=total_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))+
  labs(title="similar snps:recount hit snps (p<0.05) separated by recount total")

ggplot(ase_same %>% filter(BINOM_P<0.05 & p_val>0.05 ))+
  geom_point(aes(x= REF_RATIO, y=ref_ratio,color=total_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))+
  labs(title="similar snps:gtex hit (p<0.05) but not recount (p>0.05) hit snps separated by recount total")

ggplot(ase_same %>% filter(BINOM_P>0.05 & p_val<0.05))+
  geom_point(aes(x= REF_RATIO, y=ref_ratio,color=total_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))+
  labs(title="similar snps:not gtex hit (p>0.05) but recount (p<0.05) hit snps separated by recount total")


dev.off()

pdf(file="~/plot/ASE/ref_ratio_compare2.pdf", width = 10, height = 6)
ggplot(ase_same)+
  geom_point(aes(x= log2(REF_COUNT), y=log2(ref),color=ae_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(total_cut))+
  labs(title="similar snps:all snps separated by recount total")

ggplot(ase_same %>% filter(BINOM_P<0.05 & p_val<0.05))+
  geom_point(aes(x= log2(REF_COUNT), y=log2(ref),color=ae_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(ae_cut))+
  labs(title="similar snps:hit snps in both dataset (p<0.05) separated by recount total")

ggplot(ase_same %>% filter(BINOM_P<0.05))+
  geom_point(aes(x= log2(REF_COUNT), y=log2(ref),color=ae_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(ae_cut))+
  labs(title="similar snps:gtex hit snps (p<0.05) separated by recount total")

ggplot(ase_same %>% filter(p_val<0.05))+
  geom_point(aes(x= log2(REF_COUNT), y=log2(ref),color=ae_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(ae_cut))+
  labs(title="similar snps:recount hit snps (p<0.05) separated by recount total")

ggplot(ase_same %>% filter(BINOM_P<0.05 & p_val>0.05 ))+
  geom_point(aes(x= log2(REF_COUNT), y=log2(ref),color=ae_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(ae_cut))+
  labs(title="similar snps:gtex hit (p<0.05) but not recount (p>0.05) hit snps separated by recount total")

ggplot(ase_same %>% filter(BINOM_P>0.05 & p_val<0.05))+
  geom_point(aes(x= log2(REF_COUNT), y=log2(ref),color=ae_cut))+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  facet_wrap(vars(ae_cut))+
  labs(title="similar snps:not gtex hit (p>0.05) but recount (p<0.05) hit snps separated by recount total")


dev.off()
plot_df<-ase_same %>% select(REF_RATIO,ref_ratio,total_cut,BINOM_P,p_val) %>% pivot_longer(!c(BINOM_P,p_val,total_cut),names_to="dataset", values_to = "ratio")

pdf(file="~/plot/ASE/ref_ratio_compare_box.pdf", width = 10, height = 4)
ggplot(plot_df)+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))


ggplot(plot_df %>% filter(BINOM_P<0.05))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))

ggplot(plot_df %>% filter(p_val<0.05))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))

ggplot(plot_df %>% filter(BINOM_P<0.05 & p_val>0.05 ))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))

ggplot(plot_df %>% filter(BINOM_P>0.05 & p_val<0.05))+
  geom_boxplot(aes(x= total_cut, y=ratio,color=dataset))

dev.off()


ase_same<-ase_same %>% 
  mutate(ratio=log2(TOTAL_COUNT) - log2(total),
         mean=(log2(TOTAL_COUNT) + log2(total))/2,
         ratio_ref=log2(REF_COUNT) - log2(ref),
         mean_ref=(log2(REF_COUNT) + log2(ref))/2,
         ratio_alt=log2(ALT_COUNT) - log2(alt),
         mean_alt=(log2(ALT_COUNT) + log2(alt))/2)


pdf(file="~/plot/ASE/ae_cut_MA.pdf", width = 10, height = 4)

#ggplot(ase_same)+
# geom_point(aes(x=mean,y=ratio))
for(i in levels(ase_same$ae_cut)){
  p=ggplot(ase_same %>% filter(ae_cut==i),aes(x=mean,y=ratio))+
    geom_point(alpha=0.4)+
    geom_smooth()+
    geom_point(data=ase_same%>% filter(ae_cut==i,BINOM_P>0.01 & p_val<0.01 ),
               aes(x=mean,y=ratio), color="red", alpha=0.5)+
    labs(title=paste0("total MA plot: recount hit, no gtex hit AE=", i))
  print(p)
}
dev.off()

pdf(file="~/plot/ASE/ae_cut_MA2.pdf", width = 10, height = 4)

#ggplot(ase_same)+
# geom_point(aes(x=mean,y=ratio))
for(i in levels(ase_same$ae_cut)){
  p=ggplot(ase_same %>% filter(ae_cut==i))+
    geom_point(aes(x=mean,y=ratio), alpha=0.4)+
    geom_smooth(aes(x=mean,y=ratio))+
    geom_point(data=ase_same%>% filter(ae_cut==i,BINOM_P<0.01 & p_val>0.01 ),
               aes(x=mean,y=ratio), color="red", alpha=0.5)+
    labs(title=paste0("total MA plot:gtex hit, no recount hit AE=", i))
  print(p)
}
dev.off()

pdf(file="~/plot/ASE/ae_cut_MA_ref.pdf", width = 10, height = 4)

#ggplot(ase_same)+
 # geom_point(aes(x=mean,y=ratio))
for(i in levels(ase_same$ae_cut)){
p=ggplot(ase_same %>% filter(ae_cut==i),aes(x=mean_ref,y=ratio_ref))+
  geom_point(alpha=0.4)+
  geom_smooth()+
  geom_point(data=ase_same%>% filter(ae_cut==i,BINOM_P>0.01 & p_val<0.01 ),
             aes(x=mean_ref,y=ratio_ref), color="red", alpha=0.5)+
  labs(title=paste0("ref MA plot: recount hit, no gtex hit AE=", i))
print(p)
}
dev.off()

pdf(file="~/plot/ASE/ae_cut_MA2_ref.pdf", width = 10, height = 4)

#ggplot(ase_same)+
# geom_point(aes(x=mean,y=ratio))
for(i in levels(ase_same$ae_cut)){
  p=ggplot(ase_same %>% filter(ae_cut==i))+
    geom_point(aes(x=mean_ref,y=ratio_ref), alpha=0.4)+
    geom_smooth(aes(x=mean_ref,y=ratio_ref))+
    geom_point(data=ase_same%>% filter(ae_cut==i,BINOM_P<0.01 & p_val>0.01 ),
               aes(x=mean_ref,y=ratio_ref), color="red", alpha=0.5)+
    labs(title=paste0("ref MA plot:gtex hit, no recount hit AE=", i))
  print(p)
}
dev.off()

pdf(file="~/plot/ASE/ae_test2.pdf", width = 10, height = 4)

ggplot(ase_same, aes(x=log2(ref),y=log2(alt)))+
  geom_point(alpha=0.2)+
  geom_point(data=ase_same%>% filter(p_val<0.01 ),
             aes(x=log2(ref),y=log2(alt)), color="red", alpha=0.5)


ggplot(ase_same, aes(x=log2(REF_COUNT),y=log2(ALT_COUNT)))+
  geom_point(alpha=0.2)+
  geom_point(data=ase_same%>% filter(p_val<0.01 ),
             aes(x=log2(REF_COUNT),y=log2(ALT_COUNT)), color="red", alpha=0.5)
for(i in levels(ase_same$ae_cut)){
  p=ggplot(ase_same %>% filter(ae_cut==i),aes(x=log2(ref),y=log2(alt)))+
    geom_point(alpha=0.4)+
    geom_point(data=ase_same%>% filter(ae_cut==i,BINOM_P>0.05 & p_val<0.05 ),
               aes(x=log2(ref),y=log2(alt)), color="red", alpha=0.5)+
    labs(title=paste0("total MA plot: recount hit, no gtex hit AE=", i))
  print(p)
}
dev.off()



t<-ase_same %>%  filter(BINOM_P>0.01 & p_val<0.01 )
t<-as.data.frame(t)
median(ase_same$ref_ratio)
binom.test(63,98,0.515458)$p.value
t[2,c(1:7,12:16,19,22:23)]
0.515458
dim(t)
log2((64 + 1)/(62 + 1))
pdf(file="~/plot/ASE/outlier.pdf", width = 10, height = 4)

ggplot(ase_same, aes(x=log2(TOTAL_COUNT),y=log2(total)))+
  geom_point(alpha=0.2)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")

ggplot(ase_same, aes(x=log2(REF_COUNT),y=log2(ref)))+
  geom_point(alpha=0.2)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")

  dev.off()
pdf(file="~/plot/ASE/ae_hit.pdf", width = 10, height = 4)

ggplot(ase_same, aes(x=AE,y=ae))+
  geom_point(alpha=0.2)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  labs(title="all SNPs")
  
for(i in levels(ase_same$total_cut)){
  p=ggplot(ase_same %>% filter(total_cut==i),aes(x=AE,y=ae))+
    geom_point(alpha=0.1)+geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
    geom_point(data=ase_same%>% filter(total_cut==i,BINOM_P>0.05 & p_val<0.05 ),
               aes(x=AE,y=ae), color="red", alpha=0.3)+
    geom_point(data=ase_same%>% filter(total_cut==i,BINOM_P<0.05 & p_val>0.05 ),
               aes(x=AE,y=ae), color="blue", alpha=0.3)+
    geom_point(data=ase_same%>% filter(total_cut==i,BINOM_P<0.05 & p_val<0.05 ),
               aes(x=AE,y=ae), color="green", alpha=0.3)+
    labs(title=paste0("AE comparison: red=only hit in recount, blue=only hit in gtex, green=both"),
         subtitle= paste0("total_cut=", i))
  print(p)
}
dev.off()

pdf(file="~/plot/ASE/refRatio_hit.pdf", width = 10, height = 4)

ggplot(ase_same, aes(x=REF_RATIO,y=ref_ratio))+
  geom_point(alpha=0.2)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  labs(title="All SNPs")

for(i in levels(ase_same$total_cut)){
  p=ggplot(ase_same %>% filter(total_cut==i),aes(x=REF_RATIO,y=ref_ratio))+
    geom_point(alpha=0.1)+geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
    geom_point(data=ase_same%>% filter(total_cut==i,BINOM_P>0.05 & p_val<0.05 ),
               aes(x=REF_RATIO,y=ref_ratio), color="red", alpha=0.3)+
    geom_point(data=ase_same%>% filter(total_cut==i,BINOM_P<0.05 & p_val>0.05 ),
               aes(x=REF_RATIO,y=ref_ratio), color="blue", alpha=0.3)+
    geom_point(data=ase_same%>% filter(total_cut==i,BINOM_P<0.05 & p_val<0.05 ),
               aes(x=REF_RATIO,y=ref_ratio), color="green", alpha=0.3)+
    labs(title="ref_ratio comparison: red=only hit in recount, blue=only hit in gtex, green=both",
         subtitle= paste0("total_cut=", i))
  print(p)
}
dev.off()

pdf(file="~/plot/ASE/ref_hit.pdf", width = 10, height = 4)

ggplot(ase_same, aes(x= log2(REF_COUNT), y=log2(ref)))+
  geom_point(alpha=0.2)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  labs(title="All SNPs")

for(i in levels(ase_same$total_cut)){
  p=ggplot(ase_same %>% filter(total_cut==i),aes(x= log2(REF_COUNT), y=log2(ref)))+
    geom_point(alpha=0.1)+geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
    geom_point(data=ase_same%>% filter(total_cut==i,BINOM_P>0.05 & p_val<0.05 ),
               aes(x= log2(REF_COUNT), y=log2(ref)), color="red", alpha=0.3)+
    geom_point(data=ase_same%>% filter(total_cut==i,BINOM_P<0.05 & p_val>0.05 ),
               aes(x= log2(REF_COUNT), y=log2(ref)), color="blue", alpha=0.3)+
    geom_point(data=ase_same%>% filter(total_cut==i,BINOM_P<0.05 & p_val<0.05 ),
               aes(x= log2(REF_COUNT), y=log2(ref)), color="green", alpha=0.3)+
    labs(title="ref_ratio comparison: red=only hit in recount, blue=only hit in gtex, green=both",
         subtitle= paste0("total_cut=", i))
  print(p)
}
dev.off()

pdf(file="~/plot/ASE/ref_ratio_gtex.pdf", width = 10, height = 4)
ggplot(ase_same%>% filter(BINOM_P_ADJUSTED>0.05)) + 
  geom_boxplot(aes(x=sample_id, y=REF_RATIO),position=position_dodge(1))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+labs(title="GTEx non significant")

ggplot(ase_same %>% filter(BINOM_P_ADJUSTED<0.05)) + 
  geom_boxplot(aes(x=sample_id, y=REF_RATIO),position=position_dodge(1))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+labs(title="GTEx significant")

dev.off()

pdf(file="~/plot/ASE/ae_gtex_sig.pdf", width = 10, height = 4)
ggplot(ase_same%>% filter(BINOM_P_ADJUSTED>0.05)) + 
  geom_point(aes(x=AE, y=ae),alpha=0.5)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  labs(title="gtex none significant")

ggplot(ase_same %>% filter(BINOM_P_ADJUSTED<0.05)) + 
  geom_boxplot(aes(x=sample_id, y=REF_RATIO),alpha=0.5)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype="dotdash")+
  labs(title="gtex significant")
dev.off()

x<-ase_join%>% filter(sample_id %in% unique(ase_join$sample_id)[4:6])

pdf(file="~/plot/ASE/ASE_sign.pdf", width = 10, height = 4)
ggplot(x) + 
  geom_point(aes(x=log2(ref), y=log2(alt)),alpha=0.2)+
  geom_point(data= x%>% filter(p_val<0.05),
             aes(x=log2(ref), y=log2(alt)),alpha=0.4, color="red")+
  labs(title="Recount p-val")+
  facet_wrap(vars(sample_id))

ggplot(x) + 
  geom_point(aes(x=log2(REF_COUNT), y=log2(ALT_COUNT)),alpha=0.2)+
  geom_point(data= x%>% filter(BINOM_P<0.05),
             aes(x=log2(REF_COUNT), y=log2(ALT_COUNT)),alpha=0.4,color="red")+
  labs(title="GTEx p-val")+
  facet_wrap(vars(sample_id))

ggplot(x) + 
  geom_point(aes(x=log2(ref), y=log2(alt)),alpha=0.2)+
  geom_point(data= x%>% filter(q_val<0.05),
             aes(x=log2(ref), y=log2(alt)),alpha=0.4, color="red")+
  labs(title="Recount q-val")+
  facet_wrap(vars(sample_id))

ggplot(x) + 
  geom_point(aes(x=log2(REF_COUNT), y=log2(ALT_COUNT)),alpha=0.2)+
  geom_point(data= x%>% filter(BINOM_P_ADJUSTED<0.05),
             aes(x=log2(REF_COUNT), y=log2(ALT_COUNT)),alpha=0.4,color="red")+
  labs(title="GTEx q-val")+
  facet_wrap(vars(sample_id))

dev.off()


x%>% group_by(sample_id) %>% filter(BINOM_P<0.05) %>% summarize(n())


x%>% group_by(sample_id) %>% filter(BINOM_P_ADJUSTED<0.05) %>% summarize(n())
x%>% group_by(sample_id) %>% filter(q_val<0.05) %>% summarize(n())
