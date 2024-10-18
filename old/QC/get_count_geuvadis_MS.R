setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)

#meta_sra<-read.csv("~/ASE/data/metaSRA-runs.csv")
#sra_geno<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")
#recount3_metadata<-fread("/dcs04/hansen/data/recount_genotype/PCA/SRA/Recount3_metadata.tsv", header= T, sep = "\t",quote="")
#recount3_metadata<-recount3_metadata[,c(2:5,163)]
#id<-which(sra_geno$sample_id%in%meta_sra$sra_run_id)
geu<-read.csv("~/ASE/data/Geuvadis_metadata.csv")
geu_geno_meta<- "/dcs04/hansen/data/recount_genotype/pipeline/ERP001942/predict_genotype_accuracy/"


#Gtex paper suggests removing vriants in HLA genes:
#his BED file contains genomic positions that we have identified as either showing bias in simulations or having a UCSC mappability score < 50. 
#Variants that fall into these positions are used for phasing, but not for generating haplotypic counts to avoid problems with mapping bias.
#https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
#Downloaded blacklist of snps from gtex (https://github.com/secastel/phaser/blob/master/phaser/README.md)
bad_snp<-read.table("~/plot/ASE/hg38_haplo_count_blacklist.chr.bed", sep="\t")
bad_snp_gr<-makeGRangesFromDataFrame(bad_snp,seqnames="V1",start.field ="V2",end.field = "V3")





k=1
for(k in 1:length(geu$sample_id_rep)){
  sample_id<-geu$sample_id_rep[k]
  study<-geu$study[k]
  print(k)
  geno_dir<-paste0(geu_geno_meta, sample_id,"_predGenotypes_w_accuracy.csv.gz" )
  if(!file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase_MS.rda")))
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
      
        save(ase_df, file=paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase_MS.rda") )
      }
    }}}

load("~/test/met.rda")

#-------------------
#plot
#-------------------
xx=0
for(k in 1:nrow(geu)){
  print(k) 
  sample_id<-geu$sample_id_rep[k]
  study<-geu$study[k]
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda"))#named ase_all
    
    #load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase_MS.rda"))#named ase_df
    #all.equal(ase_all,ase_df)
    if(ase_all$num_dupl[1]>0){
      xx<-xx+1
    }
}


ase_all_liver$pred.type<-metadata_files$pred.type[match(ase_all_liver$sample_id, metadata_files$external_id)]
ase_all_liver$study<-metadata_files$study[match(ase_all_liver$sample_id, metadata_files$external_id)]

table(ase_all_liver$pred.type)

pdf(file="~/plot/ASE/SRA_normalLiver2.pdf", width = 10, height = 4)
# ggplot(ase_all_liver, aes(x=sample_id,y=ref_ratio))+
#   geom_boxplot()+
#   geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
#   geom_hline(yintercept = 0.45, color="blue",linetype="dotdash")+
#   geom_hline(yintercept = 0.55, color="purple",linetype="dotdash")+
#   labs(title="SRA normal liver samples")
# 
x<- ase_all_liver %>% group_by(study,sample_id,pred.type ) %>% summarize(is.median=median(ref_ratio)) %>% ungroup()
ggplot(x %>% filter(pred.type=="rna-seq"),aes(is.median, fill=study))+
  geom_histogram(alpha=0.5)+
  geom_vline(xintercept = 0.5, color="red",linetype="dotdash")+
  #facet_wrap(vars(pred.type))+
  theme(legend.position="none")+
  labs(title="SRA ASE analysis (2680 samples): all sample ref-ratio median density plot",
       subtitle="I chose the 'healthy/normal' liver samples from SRA and looked at the ref-ratio median of samples. Separated samples based on seq-type")

s<-sample(unique(ase_all_liver$sample_id[ase_all_liver$pred.type=="rna-seq"]),6)
dat<-ase_all_liver %>% filter(pred.type=="rna-seq",sample_id %in% s)

ggplot()+
  geom_point(data=dat %>% filter(q_val>0.05), aes(x=log2(ref), y=log2(alt)),alpha=0.3)+
  geom_point(data=dat %>% filter(q_val<0.05),aes(x=log2(ref), y=log2(alt)),alpha=0.5, color="red")+
  facet_wrap(vars(sample_id))+
  labs(title="bulk SRA ASE analysis (6 samples at random): log2(ref) vs log2(alt)",
       subtitle="red dots are q_val<0.05")

s<-sample(unique(ase_all_liver$sample_id[ase_all_liver$pred.type=="scrna-seq"]),6)
dat<-ase_all_liver %>% filter(pred.type=="scrna-seq",sample_id %in% s)
ggplot()+
  geom_point(data=dat %>% filter(q_val>0.05), aes(x=log2(ref), y=log2(alt)),alpha=0.3)+
  geom_point(data=dat %>% filter(q_val<0.05),aes(x=log2(ref), y=log2(alt)),alpha=0.5, color="red")+
  facet_wrap(vars(sample_id))+
  labs(title="scRNA-seq SRA ASE analysis (6 samples at random): log2(ref) vs log2(alt)",
       subtitle="red dots are q_val<0.05")

dev.off()

pdf(file="~/plot/ASE/SRA_normalLiver2.pdf", width = 10, height = 4)
x<-ase_all_liver%>% filter(pred.type=="rna-seq")
for(i in 1:length(unique(x$study))){
  print(i)
  p=ggplot(x%>% filter(study== unique(x$study)[i]), aes(x=sample_id,y=ref_ratio))+
    geom_boxplot()+
    geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
    geom_hline(yintercept = 0.45, color="blue",linetype="dotdash")+
    geom_hline(yintercept = 0.55, color="purple",linetype="dotdash")+
    labs(title="SRA normal liver samples in one study: bulk RNA-seq",
         subtitle=paste0("study:",unique(x$study)[i]))
  print(p)
}
dev.off()

pdf(file="~/plot/ASE/SRA_normalLiver_sig.pdf", width = 10, height = 4)

s<-sample(unique(ase_all_liver$sample_id[ase_all_liver$pred.type=="rna-seq"]),6)
dat<-ase_all_liver %>% filter(pred.type=="rna-seq",study == "SRP095921") #SRP217753

for(i in 1:length(unique(dat$sample_id))){
  p=ggplot()+
    geom_point(data=dat %>% filter(sample_id==unique(dat$sample_id)[i],q_val>0.05), aes(x=log2(ref), y=log2(alt)),alpha=0.3)+
    geom_point(data=dat %>% filter(sample_id==unique(dat$sample_id)[i],q_val<0.05),aes(x=log2(ref), y=log2(alt)),alpha=0.5, color="red")+
    labs(title="bulk SRA ASE analysis (6 samples at random): log2(ref) vs log2(alt)",
         subtitle=paste0("study=",unique(dat$sample_id)[i],", red dots are q_val<0.05"))
  print(p)
}
dev.off()

pdf(file="~/plot/ASE/mono_allelic_filter_geu.pdf", width = 10, height = 4)
for(i in unique(dat$sample_id)[1:3]){
  print(i)
  x<-dat %>% filter(sample_id==i)%>% rowwise() %>% mutate(mine_allele=min(alt,ref))
  p=ggplot(x)+
    geom_point(aes(x=total, y=mine_allele))+xlim(0,100)+ylim(0,100)
  print(p)
}
dev.off()
