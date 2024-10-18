
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


sra<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/metadata/all_SRA.csv")


for(k in 1:150){
  print(k)
  sample_id<-sra$sample_id[k]
  study<-sra$study[k]
  if(file.exists(file=paste0("~/hansen_lab/ASE/test_ASE/", sample_id, ".rds") ))
  {
  if(!file.exists(file=paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda") ))
  {
    print(sample_id)
    ase_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, ".rds"))
    ase_df$ref_ratio<- ase_df$ref/ase_df$total
    ase_df$aFC<-log2((ase_df$alt+1)/(ase_df$ref+1))
    
    #Make Granges:
    ase_filt_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="pos",end.field = "pos", keep.extra.columns = T)
      
      #Remove blacklist snps from our granges:
      ov<-findOverlaps(ase_filt_gr,bad_snp_gr)
      print(paste0("number of blacklist: ", length(ase_filt_gr)-length(ase_filt_gr[-unique(queryHits(ov))])))
      ase_filt_gr<-ase_filt_gr[-unique(queryHits(ov))]
      
      #obtain the median of the ref ratio to be used as the null p value
      is.median<-median(ase_filt_gr$ref_ratio)
      ase_filt_gr$p_val = apply(as.data.frame(ase_filt_gr)[,c("ref","alt")], 1, function(x) binom.test(x[1],(x[1]+x[2]),p=is.median)$p.value)
      # perform multiple testing correction with FDR
      ase_filt_gr$q_val = p.adjust(ase_filt_gr$p_val, method = "fdr")
      
      
        ase_all<- as_tibble(ase_filt_gr)
      
    #order the df and save the data:
    colnames(ase_all)[1:2]<-c("chr","pos")
    ase_all[,3:5]<-NULL
    
    save(ase_all, file=paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda") )
  }}}



