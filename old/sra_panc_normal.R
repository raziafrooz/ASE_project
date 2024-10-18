setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)

meta_sra<-read.csv("~/ASE/data/metaSRA-runs-panc.csv")
sra_geno<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")
recount3_metadata<-fread("/dcs04/hansen/data/recount_genotype/PCA/SRA/Recount3_metadata.tsv", header= T, sep = "\t",quote="")
recount3_metadata<-recount3_metadata[,c(2:5,163)]
sra_geno$seq_type<-recount3_metadata$seq_type[match(sra_geno$sample_id,recount3_metadata$external_id)]
sra_geno<-sra_geno[which(sra_geno$sample_id%in%meta_sra$sra_run_id),]
sra_geno<-sra_geno[which(sra_geno$seq_type=="bulk"),]

#Gtex paper suggests removing vriants in HLA genes:
#his BED file contains genomic positions that we have identified as either showing bias in simulations or having a UCSC mappability score < 50. 
#Variants that fall into these positions are used for phasing, but not for generating haplotypic counts to avoid problems with mapping bias.
#https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
#Downloaded blacklist of snps from gtex (https://github.com/secastel/phaser/blob/master/phaser/README.md)
bad_snp<-read.table("~/plot/ASE/hg38_haplo_count_blacklist.chr.bed", sep="\t")
bad_snp_gr<-makeGRangesFromDataFrame(bad_snp,seqnames="V1",start.field ="V2",end.field = "V3")


> unique(sra_geno$study)
[1] "ERP003613" "ERP113252" "SRP013565" "SRP014739" "SRP015640" "SRP039090"
[7] "SRP042282" "SRP048556" "SRP052057" "SRP055513" "SRP077921" "SRP219514"



for(k in 1:nrow(sra_geno)){
  print(k)
  sample_id<-sra_geno$sample_id[k]
  study<-sra_geno$study[k]
  geno_dir<-sra_geno$genotypedSamples[k]
  if(!file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda")))
  {
    
    ase_df<-as_tibble(read.csv(geno_dir))%>% 
      filter(pred_genotype==2, coverage>=8)
    
    if(nrow(ase_df)>10){
      ase_df <-ase_df %>%
        mutate(alt=(sqrt(2^((2*S)-M))-1),ref=((2^M)*(alt+1))-1) %>% 
        select(!c(M,S,pred_genotype,predicted_accuracy))
      colnames(ase_df)[c(2,4)]<-c("pos","total")
      #ase_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, ".rds"))
      ase_df$ref_ratio<- ase_df$ref/ase_df$total
      ase_df$aFC<-log2((ase_df$alt+1)/(ase_df$ref+1))
      
      #Make Granges:
      ase_filt_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="pos",end.field = "pos")
      
      #Remove blacklist snps from our granges:
      ov<-findOverlaps(ase_filt_gr,bad_snp_gr)
      print(paste0("number of blacklist: ", length(unique(queryHits(ov)))))
      if(length(unique(queryHits(ov)))>0){
        ase_df<-ase_df[-unique(queryHits(ov)),]
        
        #obtain the median of the ref ratio to be used as the null p value
        is.median<-median(ase_df$ref_ratio)
        print(paste0("median is ",is.median))
        ase_df$p_val = apply(ase_df[,c("ref","alt")], 1, function(x) binom.test(x[1],(x[1]+x[2]),p=is.median)$p.value)
        # perform multiple testing correction with FDR
        ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")
        
        save(ase_df, file=paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda") )
      }
    }}}

load("~/test/met.rda")

#-------------------
#plot
#-------------------

for(k in 1:nrow(sra_geno)){
  print(k)
  sample_id<-sra_geno$sample_id[k]
  study<-sra_geno$study[k]
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda")))
  {
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda"))
    ase_df$sample_id<-sample_id
    ase_df$study<-study
    if(k==1){
      ase_all_panc<-ase_df
    }
    ase_all_panc<-rbind(ase_all_panc,ase_df)
  }}

unique(ase_all_panc$sample_id)
pdf(file="~/plot/ASE/SRA_normalPanc_all.pdf", width = 10, height = 4)

x<- ase_all_panc %>% group_by(study,sample_id ) %>% summarize(is.mean=mean(ref_ratio)) %>% ungroup()
ggplot(x,aes(y=is.mean, x=study))+
  geom_point(alpha=0.5)+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  #facet_wrap(vars(pred.type))+
  theme(legend.position="none")+
  labs(title="SRA healthy pacreatic ASE analysis (133 samples): all sample ref-ratio mean",
       subtitle="Each dot shows a sample. Study is on x axis. Only bulk RNA-seq")

dev.off()

pdf(file="~/plot/ASE/SRA_normalPanc.pdf", width = 10, height = 4)
x<-ase_all_panc
for(i in 1:length(unique(x$study))){
  print(i)
  p=ggplot(x%>% filter(study== unique(x$study)[i]), aes(x=sample_id,y=ref_ratio))+
    geom_boxplot()+
    geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
    geom_hline(yintercept = 0.45, color="blue",linetype="dotdash")+
    geom_hline(yintercept = 0.55, color="purple",linetype="dotdash")+
    labs(title="SRA normal Pancreatic samples in one study: bulk RNA-seq",
         subtitle=paste0("study:",unique(x$study)[i]))
  print(p)
}
dev.off()

s<-sample(unique(ase_all_panc$sample_id),6)
dat<-ase_all_panc %>% filter(sample_id %in% s) #SRP217753

pdf(file="~/plot/ASE/SRA_normalPanc_sig.pdf", width = 10, height = 4)

for(i in 1:length(unique(dat$sample_id))){
  p=ggplot()+
    geom_point(data=dat %>% filter(sample_id==unique(dat$sample_id)[i],q_val>0.05), aes(x=log2(ref), y=log2(alt)),alpha=0.3)+
    geom_point(data=dat %>% filter(sample_id==unique(dat$sample_id)[i],q_val<0.05),aes(x=log2(ref), y=log2(alt)),alpha=0.5, color="red")+
    labs(title="bulk SRA ASE analysis (6 samples at random): log2(ref) vs log2(alt)",
         subtitle=paste0("study=",unique(dat$sample_id)[i],", red dots are q_val<0.05"))
  print(p)
}
dev.off()

pdf(file="~/plot/ASE/mono_allelic_filter_panc.pdf", width = 10, height = 4)
for(i in unique(ase_all_panc$sample_id)[1:10]){
  print(i)
  x<-ase_all_panc %>% filter(sample_id==i)%>% rowwise() %>% mutate(mine_allele=min(alt,ref))
  p=ggplot(x)+
    geom_point(aes(x=total, y=mine_allele))+xlim(0,100)+ylim(0,100)+
    labs(title="First 10 Pancreatic samples: mono-allele check")
  print(p)
}
dev.off()
