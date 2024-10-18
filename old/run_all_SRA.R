setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)

sra_geno<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")
recount3_metadata<-fread("/dcs04/hansen/data/recount_genotype/PCA/SRA/Recount3_metadata.tsv", header= T, sep = "\t",quote="")
recount3_metadata<-recount3_metadata[,c(2:5,163)]
sra_geno$seq_type<-recount3_metadata$seq_type[match(sra_geno$sample_id,recount3_metadata$external_id)]
sra_geno<-sra_geno[which(sra_geno$seq_type=="bulk"),]

#Gtex paper suggests removing vriants in HLA genes:
#his BED file contains genomic positions that we have identified as either showing bias in simulations or having a UCSC mappability score < 50. 
#Variants that fall into these positions are used for phasing, but not for generating haplotypic counts to avoid problems with mapping bias.
#https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
#Downloaded blacklist of snps from gtex (https://github.com/secastel/phaser/blob/master/phaser/README.md)
bad_snp<-read.table("~/plot/ASE/hg38_haplo_count_blacklist.chr.bed", sep="\t")
bad_snp_gr<-makeGRangesFromDataFrame(bad_snp,seqnames="V1",start.field ="V2",end.field = "V3")



number_of_bad_ref<-0
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
        dplyr::select(!c(M,S,pred_genotype,predicted_accuracy))
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
        if(is.median<0.4 | is.median>0.6){
          number_of_bad_ref<-number_of_bad_ref+1
        }
        print(paste0("median is ",is.median))
        ase_df$p_val = apply(ase_df[,c("ref","alt")], 1, function(x) {
          binom.test(round(x[1],1),round((x[1]+x[2]),1),p=is.median)$p.value})
        # perform multiple testing correction with FDR
        ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")
        
        save(ase_df, file=paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda") )
      }
    }}
  print(paste0(number_of_bad_ref, " samples have bad ref_ratio"))
}


