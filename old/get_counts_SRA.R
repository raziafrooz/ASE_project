setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

sra<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/metadata/all_SRA.csv")
sra[1,]
recount3_chr_mapping <- "/dcl02/lieber/ajaffe/recount-pump/recount3.alts.chromosome_mappings.tsv"
snp_gr<-readRDS("/dcs04/hansen/data/recount_genotype/biallelic_snps/biallelic_SNP_gr.rds")
sra_geno<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")

k=1
#########################
#SRA samples:
#########################

coverage_cutoff=8
for(k in 1:500){
  print(k)
  sample_id<-sra$sample_id[k]
  study<-sra$study[k]
  if(!file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, ".rds")))
  {
  geno<-as_tibble(read.csv(sra_geno$genotypedSamples[k]))
  geno <-geno %>% filter(pred_genotype==2, coverage>=8)
  if(nrow(geno)>0){
  
  het_gr<-makeGRangesFromDataFrame(geno,seqnames="chr",start.field ="start",end.field = "start")

  
  
    print(sample_id)
    bigWig_path <- sra$total[k]
    cat("Loading in: bigWig", "\n")
    
    bigwig <- tryCatch(
      {
        import(bigWig_path, format = "bigwig")
      },
      error=function(cond) {
        message(paste("Error loading bigwig file: ", bigWig_path))
        message(cond)
        return(NA)
      },
      warning=function(cond) {
        message(paste("Warning loading bigwig file: ", bigWig_path))
        message(cond)
        return(NA)
      },
      finally={}
    )    
    if(all(is.na(bigwig))) {
      return(NA)
    }
    overlap_loci <- findOverlaps(snp_gr, bigwig)
    bigwig_count <- bigwig$score[subjectHits(overlap_loci)]
    filter_idx <- which(bigwig_count >= coverage_cutoff)
    #bigwig_count <- bigwig_count[filter_idx]
    filtered_snps_gr <- snp_gr[queryHits(overlap_loci)[filter_idx]]
    filtered_snps_gr$total<-bigwig_count[filter_idx]
    #Remove homozygous snps from our granges:
    ov<-findOverlaps(filtered_snps_gr,het_gr)
    
    filtered_snps_gr<-filtered_snps_gr[queryHits(ov)]
    #geno_data<-readRDS(geno_met$allGenotypesOutput[geno_met$study==study])
    #for(i in 1:length(unique(geno_data$sample_id_rep))){
     # print(i)
      #sam<-unique(geno_data$sample_id_rep)[i]
      ##start the analysis with 1 indv, only select hetrozyogus locations, total read cound > 8)
      #geno<-geno_data %>% filter(sample_id_rep==sam, pred_genotype==2, coverage>=8) %>% select(c(chr,start,AF,coverage,pred_genotype, true_genotype, sample_id_rep,predicted.values.prob ))# SNPs= 23 328
      #geno_gr<-makeGRangesFromDataFrame(geno, start.field = "start", end.field = "start", keep.extra.columns = T)
      #genome(geno_gr)<-"hg38"
      
      
      #############
      #Read in the alt file
      ##############
    
      temp="~/test"
      temp_altFile <- paste0(temp, sample_id, ".alt.temp.csv")
      system(paste0("zstdcat ", sra$alt[k], " > ", temp_altFile))
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
      
      ov <- findOverlaps(alt_count_gr,filtered_snps_gr)
      
      #################################################
      #this is not exclusively biallelic!!!!!!!!!!!! come back to it and filter
      #################################################
      #Subset `alt_count_gr`, `alt_count` to the positions from SNPs, and compute 
      #its own unique key.
      alt_count_gr <- alt_count_gr[queryHits(ov)]
      alt_count <- alt_count[queryHits(ov)]
      alt_key <- paste(as.character(seqnames(alt_count_gr)), start(alt_count_gr), 
                       filtered_snps_gr$ref_seq[subjectHits(ov)], alt_count$alt, sep = "_")
      #Even though `alt_key` is in the same positions as `filtered_snps_key`, many of the 
      #`alt_key` entries refer to alternate alleles that we are not tracking. 
      #We need to match a second time, this time using the entire key.

      filtered_snps_key <- paste(as.character(seqnames(filtered_snps_gr)), start(filtered_snps_gr), 
                                 filtered_snps_gr$ref_seq, filtered_snps_gr$alt_seq, sep = "_")
      idx <- match(filtered_snps_key, alt_key)
      snps_key_idx <- which(!is.na(idx))
      alt_key_idx <- idx[!is.na(idx)]
      #Construct final `alt_count` relative to `filtered_snps_gr`
      #final_alt_count <- rep(0, length(filtered_snps_gr))
      filtered_snps_gr$alt<-0
      filtered_snps_gr$alt[snps_key_idx] <- alt_count[alt_key_idx]$N
      
      
      #filtered_snps_gr<-filtered_snps_gr[!is.na(filtered_snps_gr$alt)]
      #We sometimes have cases where the alt counts > coverage counts (bigWig): 
      #the alt reads were processed to keep overlapping pair-end reads
      #whereas the coverage counts (bigWig) did not keep overlapping pair-end reads. 
      #this is an ad hoc way to deal with negative counts in `ref_mtx`.
      print(paste0("duplicated count is:", sum(filtered_snps_gr$total<filtered_snps_gr$alt)))
      filtered_snps_gr[which(filtered_snps_gr$alt > filtered_snps_gr$total),] <-NULL ###for now just exclude them
      
      #Get ref count:
      filtered_snps_gr$ref<-filtered_snps_gr$total-filtered_snps_gr$alt
      
      
      ase_df<-tibble(chr=as.character(seqnames(filtered_snps_gr)),pos=start(filtered_snps_gr),
                      sample_id=sample_id,
                      study=study,
                      ref=filtered_snps_gr$ref, 
                      alt= filtered_snps_gr$alt, 
                      total=filtered_snps_gr$total)
      
      IS.MEDIAN<-median(ase_df$ref/ase_df$total)
      print(paste0("Median is ", IS.MEDIAN ))
      warning(if(!0.4<IS.MEDIAN &IS.MEDIAN <0.6){print("Median is not 0.5")}) 
      
      saveRDS(ase_df,file = paste0("~/hansen_lab/ASE/test_ASE/", sample_id, ".rds"))

    }
}}
  
  


