library(optparse)

option_list <- list(
  make_option(c("--geno_dir"), type = "character",
              help = "Input: String: path to the genotyped file for a sample."),
  make_option(c("--result"), type = "character",
              help = "Input: String: Path to final ASE result."))

opt <- parse_args(OptionParser(option_list = option_list))

if (length(opt$geno_dir) == 0 |
    length(opt$result) == 0 ) {
  stop("Not all arguments provided. Check --help for description.")
}



setwd("~/ASE/")
suppressPackageStartupMessages({
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)
})




  cat("loading in the genotype file \n")
# Read in the data and only select heterozygous SNPs with coverage 8 or higher
    ase_df<-as_tibble(read.csv(opt$geno_dir))%>% 
      filter(pred_genotype==2, coverage>=8) %>% select(-M,-S)
    
    if(nrow(ase_df)==0){
      print(paste0("sample has ", nrow(ase_df), "heterozygous SNPS: remove"))
      fwrite(ase_df, file=opt$result)
      q()
    }else{
#------------------
#read in files:
      #Gtex paper suggests removing vriants in HLA genes:
      #his BED file contains genomic positions that we have identified as either showing bias in simulations or having a UCSC mappability score < 50. 
      #Variants that fall into these positions are used for phasing, but not for generating haplotypic counts to avoid problems with mapping bias.
      #https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
      #Downloaded blacklist of snps from gtex (https://github.com/secastel/phaser/blob/master/phaser/README.md)
      HLA_snp<-read.table("/users/arazi/hansen_lab/dwl_files/ASE_filter/hg38_haplo_count_blacklist.chr.bed", sep="\t")
      HLA_snp_gr<-makeGRangesFromDataFrame(HLA_snp,seqnames="V1",start.field ="V2",end.field = "V3")
      
      #I found a bed file including the problematic sites from ENCODE
      #https://github.com/Boyle-Lab/Blacklist/tree/master
      black_list<-import("/users/arazi/hansen_lab/dwl_files/ASE_filter/hg38-blacklist.v2.bed.gz") 
      
      #https://bismap.hoffmanlab.org
      mappability<-import("/users/arazi/hansen_lab/dwl_files/ASE_filter/k50.Unique.Mappability.bb")
#-----------------      
      
      ase_df$ref_ratio<- ase_df$ref_count/ase_df$coverage
      ase_df$alt_count<-ase_df$coverage - ase_df$ref_count
      #Make Granges:
      ase_filt_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")
      
      #Remove blacklist snps from our granges:
      ov<-findOverlaps(ase_filt_gr,HLA_snp_gr)
      ase_df<-ase_df[-unique(queryHits(ov)),]
      ase_filt_gr<-ase_filt_gr[-unique(queryHits(ov)),]
      
      ov<-findOverlaps(ase_filt_gr,black_list)
      ase_df<-ase_df[-unique(queryHits(ov)),]
      ase_filt_gr<-ase_filt_gr[-unique(queryHits(ov)),]
      
      ov<-findOverlaps(ase_filt_gr,mappability)
      ase_df<-ase_df[unique(queryHits(ov)),]
      ase_filt_gr<-ase_filt_gr[unique(queryHits(ov)),]
      
      if(nrow(ase_df)==0){
        print(paste0("sample has ", nrow(ase_df), "heterozygous SNPS: remove"))
        fwrite(ase_df, file=opt$result)
        q()
      }else{
        #calculate median ref-ratio to be used as the binomial distribution p value
        is.median<-median(ase_df$ref_ratio)
        print(paste0("median is ",is.median))
        # Get the binomial p-value
        ase_df$p_val = apply(ase_df[,c("ref_count","alt_count")], 1, function(x) {
          binom.test(round(x[1],1),round((x[1]+x[2]),1),p=is.median)$p.value})
        
        # perform multiple testing correction with FDR
        ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")
        
        fwrite(ase_df, file=opt$result )
      }
    }
    