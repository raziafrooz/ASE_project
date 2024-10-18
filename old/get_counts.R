# In this analysis I am using the ASE data from GTEx
#This data was downloaded by Nick during rotations from dbGap: V8 not WASP corrected :location /dcl01/hansen/data/gtex_ase
#The GTEx paper recommends using the ase WASP corrected (allelic mapping error). Afrooz downloaded this on 08/09/2023 and is here /dcl01/hansen/data/arazi/ASE/dbGap

###For now start the analysis with ase not WASP corrected:

setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)


recount3_chr_mapping <- "/dcl02/lieber/ajaffe/recount-pump/recount3.alts.chromosome_mappings.tsv"
snp_gr<-readRDS("/dcs04/hansen/data/recount_genotype/biallelic_snps/biallelic_SNP_gr.rds")

met<-read.csv("data/GTEx_metadata.csv")
geno_met<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")

#########################
#start with kidney medulla since it only has 4 samples
#########################
tissue<-c("Liver", "Lung", "Stomach", "Pancreas")
tissue_ab<-c("LIVER", "LUNG","STMACH","PNCREAS")

for(k in 1:length(tissue)){
  tissue_name<-tissue[k]
  if(!file.exists(paste0("~/hansen_lab/ASE/test_ASE/", tissue_name, ".rds")))
    {
    print(tissue_name)
geno_data<-readRDS(geno_met$allGenotypesOutput[geno_met$study==tissue_name])
for(i in 1:length(unique(geno_data$sample_id_rep))){
  print(i)
  sam<-unique(geno_data$sample_id_rep)[i]
  ##start the analysis with 1 indv, only select hetrozyogus locations, total read cound > 8)
  geno<-geno_data %>% filter(sample_id_rep==sam, pred_genotype==2, coverage>=8) %>% select(c(chr,start,AF,coverage,pred_genotype, true_genotype, sample_id_rep,predicted.values.prob ))# SNPs= 23 328
  geno_gr<-makeGRangesFromDataFrame(geno, start.field = "start", end.field = "start", keep.extra.columns = T)
  genome(geno_gr)<-"hg38"
  
  
  #############
  #Read in the alt file
  ##############
  temp="~/test"
  temp_altFile <- paste0(temp, sam, ".alt.temp.csv")
  system(paste0("zstdcat ", met$alt[met$sample_id_rep== sam ], " > ", temp_altFile))
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
  rm(final_alt_count)
  #We sometimes have cases where the alt counts > coverage counts (bigWig): 
  #the alt reads were processed to keep overlapping pair-end reads
  #whereas the coverage counts (bigWig) did not keep overlapping pair-end reads. 
  #this is an ad hoc way to deal with negative counts in `ref_mtx`.
  
  geno_gr[geno_gr$alt_count > geno_gr$coverage] <-NULL ###for now just exclude them
  
  #Get ref count:
  geno_gr$ref_count<-geno_gr$coverage-geno_gr$alt_count
  
  
  ase_df2<-tibble(chr=as.character(seqnames(geno_gr)),pos=start(geno_gr),
                  sample_id=sam,
                  ref=geno_gr$ref_count, 
                  alt= geno_gr$alt_count, 
                  total=geno_gr$coverage,
                  pred_genotype=geno_gr$pred_genotype,
                  true_genotype=geno_gr$true_genotype)
  
  IS.MEDIAN<-median(ase_df2$ref/ase_df2$total)
  print(paste0("Median is ", IS.MEDIAN ))
  warning(if(!0.4<IS.MEDIAN &IS.MEDIAN <0.6){print("Median is not 0.5")}) 
  
  
   if(i==1){ase_df<-ase_df2} else {
    ase_df<-rbind(ase_df,ase_df2)}
}
saveRDS(ase_df,file = paste0("~/hansen_lab/ASE/test_ASE/", tissue_name, ".rds"))
  }
  

#------------------------------------------------------
#Get the true ASE hits from gtex 
#------------------------------------------------------
  if(!file.exists(paste0("~/hansen_lab/ASE/test_ASE/true_gtex_", tissue_name, ".rds")))
  {
    print(tissue_name)
    
for(i in 1:length(unique(geno_data$sample_id_rep))){
  print(i)
  sam<-unique(geno_data$sample_id_rep)[i]
  indv<-met$individual_id[met$sample_id_rep== sam ]
  gtex<-"/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_subject/"
  gtex_indv<- fread(paste0(gtex,indv,".v8.ase_table.tsv.gz" ))
  gtex_indv<-gtex_indv[gtex_indv$TISSUE_ID==tissue_ab[k],]
  gtex_indv<-as_tibble(gtex_indv[,c(1,2,9,10,11,12,15,16,20,21)])
  gtex_indv$sample_id<-sam
  if(i==1){
    true_gtex<-gtex_indv} else {
      true_gtex<-rbind(true_gtex,gtex_indv)
    }
}
saveRDS(true_gtex,file = paste0("~/hansen_lab/ASE/test_ASE/true_gtex_", tissue_name, ".rds"))
}}



