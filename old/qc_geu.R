setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)
library(ggplot2)
library(cowplot)
recount3_chr_mapping <- "/dcl02/lieber/ajaffe/recount-pump/recount3.alts.chromosome_mappings.tsv"
snp_gr<-readRDS("/dcs04/hansen/data/recount_genotype/biallelic_snps/biallelic_SNP_gr.rds")

#------------------------------------------------------
##Start with Geuvadis samples:
#------------------------------------------------------

geu<-read.csv("~/ASE/data/Geuvadis_metadata.csv")
geu_geno_meta<- "/dcs04/hansen/data/recount_genotype/pipeline/ERP001942/predict_genotype_accuracy/"

for(i in 1:length(geu$sample_id_rep)){
  print(i)
  geu_geno<- read.csv(paste0(geu_geno_meta, geu$sample_id_rep[i],"_predGenotypes_w_accuracy.csv.gz" ))
  geu_geno<-geu_geno %>% filter(pred_genotype==2) %>% select(c(chr,start,AF,coverage,pred_genotype, predicted_accuracy ))# SNPs= 9642
  geno_gr<-makeGRangesFromDataFrame(geu_geno, start.field = "start", end.field = "start", keep.extra.columns = T)
  genome(geno_gr)<-"hg38"
  
  
  #ov <- findOverlaps(snp_gr, geno_gr)
  #geno_gr <- geno_gr[subjectHits(ov)]
  
  
  
  
  #############
  #Read in the alt file
  ##############
  temp="~/test"
  temp_altFile <- paste0(temp, geu$sample_id_rep[i], ".alt.temp.csv")
  system(paste0("zstdcat ", geu$alt[i], " > ", temp_altFile))
  alt <- fread(temp_altFile, select = c(1, 2, 4))
  system(paste0("rm ", temp_altFile))
  if(nrow(alt) == 0) {
    return(NA)
  }
  colnames(alt) <- c("chr", "pos", "alt")
  #collapse alt counts. 
  alt_count <- alt[, .N, by = .(chr, pos, alt)] 
  rm(alt)
  alt_count <- alt_count[!is.na(alt_count$alt) ,]
  #`alt_count` is 0-based, so we add 1 to positions. 
  alt_count$pos <- alt_count$pos + 1
  
  #use chromosome names that are used in recount3.
  chr_mapping <- fread(recount3_chr_mapping, header = FALSE)
  alt_count$chr <- chr_mapping$V2[match(alt_count$chr, chr_mapping$V1)] 
  alt_count_gr <- GRanges(seqnames = alt_count$chr,
                          ranges = IRanges(alt_count$pos,
                                           alt_count$pos))
  
  ov <- findOverlaps(alt_count_gr,snp_gr)
  
  #################################################
  #this is not exclusively biallelic!!!!!!!!!!!! come back to it and filter
  #################################################
  #Subset `alt_count_gr`, `alt_count` to the positions from SNPs, and compute 
  #its own unique key.
  alt_count_gr <- alt_count_gr[queryHits(ov)]
  alt_count <- alt_count[queryHits(ov)]
  alt_key <- paste(as.character(seqnames(alt_count_gr)), start(alt_count_gr), 
                   snp_gr$ref_seq[subjectHits(ov)], alt_count$alt, sep = "_")
  #Even though `alt_key` is in the same positions as `filtered_snps_key`, many of the 
  #`alt_key` entries refer to alternate alleles that we are not tracking. 
  #We need to match a second time, this time using the entire key.
  ov <- findOverlaps(geno_gr,snp_gr)
  filtered_snps_key <- paste(as.character(seqnames(geno_gr)), start(geno_gr), 
                             snp_gr$ref_seq[subjectHits(ov)], snp_gr$alt_seq[subjectHits(ov)], sep = "_")
  idx <- match(filtered_snps_key, alt_key)
  snps_key_idx <- which(!is.na(idx))
  alt_key_idx <- idx[!is.na(idx)]
  #Construct final `alt_count` relative to `filtered_snps_gr`
  final_alt_count <- rep(0, length(geno_gr))
  final_alt_count[snps_key_idx] <- alt_count[alt_key_idx]$N
  
  
  geno_gr$alt_count<-final_alt_count
  #We sometimes have cases where the alt counts > coverage counts (bigWig): 
  #the alt reads were processed to keep overlapping pair-end reads
  #whereas the coverage counts (bigWig) did not keep overlapping pair-end reads. 
  #this is an ad hoc way to deal with negative counts in `ref_mtx`.
  
  geno_gr[geno_gr$alt_count > geno_gr$coverage] <-NULL ###for now just exclude them
  
  #Get ref count:
  geno_gr$ref_count<-geno_gr$coverage-geno_gr$alt_count
  #geno_gr$ref_count[geno_gr$ref_count < 0] <- 0 
  
  ase_df<-tibble(chr=as.character(seqnames(geno_gr)),pos=start(geno_gr),ref=geno_gr$ref_count, alt= geno_gr$alt_count, total=geno_gr$coverage,
                 ref_ratio=abs(geno_gr$ref_count/geno_gr$coverage),
                 p_val=NA,
                 q_val=NA)
  
  
  ########AE should be calculated this way: when you come back double check:
  AE=abs(0.5-(geno_gr$ref_count/geno_gr$coverage))
  
  for (j in 1:nrow(ase_df)){
    ase_df$p_val[j]<-binom.test(round(ase_df$ref_ratio[j] *  ase_df$total[j]),ase_df$total[j], 0.5)$p.value
  }
  #correct for multiple testing by FDR
  
  ase_df$q_val<-qvalue(ase_df$p_val,fdr.level = 0.1)$qvalues
  saveRDS(ase_df,file = paste0("~/hansen_lab/ASE/test_ASE/", geu$sample_id_rep[i], ".rds"))
}


pdf(file="~/plot/ref_ratio.pdf", width = 6, height = 4)
for(i in 1:20){
  print(i)
  ase_df<- readRDS(paste0("~/hansen_lab/ASE/test_ASE/", geu$sample_id_rep[i], ".rds"))
  ase_df$q_color<- NA
  ase_df$q_color[ase_df$q_val<0.05]<-"hit"
  
  p1<-ggplot(ase_df, aes(x=geu$sample_id_rep[i], y=ref_ratio)) + 
    geom_boxplot()+ geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+ labs(main=geu$sample_id_rep[i])
  
  p2<-ggplot(ase_df %>% mutate(log_ref=log(ref), log_alt=log(alt)), aes(x=log_ref, y=log_alt, color=q_color)) + 
    geom_point(alpha=0.4)+ geom_abline(intercept = 0, slope = 1,color="red",linetype="dotdash")+ labs(color=geu$sample_id_rep[i])
  #p3<-ggplot(ase_df) + 
  #geom_histogram(aes(x=p_val), color="black", fill=alpha("lightblue",0.4))+
  #geom_histogram(aes(x=q_val), color="black", fill=alpha("salmon",0.4))
  
  print(p1)
  print(p2)
}
dev.off()
