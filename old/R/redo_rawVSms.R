setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

met_gtex<-read.csv("data/GTEx_metadata.csv")
recount3_chr_mapping <- "/dcl02/lieber/ajaffe/recount-pump/recount3.alts.chromosome_mappings.tsv"
snps_gr<-readRDS("/dcs04/hansen/data/recount_genotype/biallelic_snps/biallelic_SNP_gr.rds")

sample_id<-"GTEX-131YS-0011-R10b-SM-5EQ5N.1"
tissue<-"Lung"
join_ase<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rds"))


alt_path<-met_gtex$alt[which(met_gtex$sample_id_rep==sample_id & met_gtex$tissue == tissue )]
bigWig_path<-met_gtex$total[which(met_gtex$sample_id_rep==sample_id & met_gtex$tissue == tissue)]

geno<-join_ase
geno_gr<-makeGRangesFromDataFrame(geno, start.field = "pos", end.field = "pos", keep.extra.columns = T)
genome(geno_gr)<-"hg38"
ov<-findOverlaps(geno_gr,snps_gr )
snp_gr<-snps_gr[subjectHits(ov)]




#----------------------------
#read alt


cat("Loading in: ", alt_path, "\n")
temp_folder="~/test"
temp_altFile <- paste0(temp_folder, sample_id, ".alt.temp.csv")
tryCatch({
  system(paste0("zstdcat ", alt_path, " > ", temp_altFile))
},
error = function(e){
  message("Error reading alt file:\n", e)
  message("'zstdcat' software is not installed")
},
warning = function(w){
  message("Warning reading alt file:\n", w)
  message("'zstdcat' software is not installed")
},
finally = {
})


alt <- data.table::fread(temp_altFile, select = c(1, 2, 4))
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
alt_count_gr <- GenomicRanges::GRanges(seqnames = alt_count$chr,
                                       ranges = IRanges(alt_count$pos,
                                                        alt_count$pos))
ov <- GenomicRanges::findOverlaps(alt_count_gr, snp_gr)
#Subset `alt_count_gr`, `alt_count` to the positions from SNPs, and compute
#its own unique key.
alt_count_gr <- alt_count_gr[queryHits(ov)]
alt_count <- alt_count[queryHits(ov)]
alt_key <- paste(as.character(seqnames(alt_count_gr)), start(alt_count_gr),
                 snp_gr$ref_seq[subjectHits(ov)], alt_count$alt, sep = "_")
#Even though `alt_key` is in the same positions as `filtered_snps_key`, many of the
#`alt_key` entries refer to alternate alleles that we are not tracking.
#We need to match a second time, this time using the entire key.
filtered_snps_key <- paste(as.character(seqnames(snp_gr)), start(snp_gr),
                           snp_gr$ref_seq, snp_gr$alt_seq, sep = "_")
idx <- match(filtered_snps_key, alt_key)
snps_key_idx <- which(!is.na(idx))
alt_key_idx <- idx[!is.na(idx)]

final_alt_count <- rep(0, length(snp_gr))
final_alt_count[snps_key_idx] <- alt_count[alt_key_idx]$N

#Compute M and S values.
ref_count <- bigwig_count - final_alt_count
#We sometimes have
all.equal(geno$alt, final_alt_count)




ov <- findOverlaps(alt_count_gr,snp_gr)
dd<-ov[which(duplicated(subjectHits(ov)))]
dd_2<-dd[which(duplicated(subjectHits(dd)))]
dd_3<-dd_2[which(duplicated(subjectHits(dd_2)))]

error1<-dd[which(!duplicated(subjectHits(dd)))]
error2<-dd_2[which(!duplicated(subjectHits(dd_2)))]
error3<-dd_3[which(!duplicated(subjectHits(dd_3)))]

geno$alt_err1<-NA
geno$alt_err1_count<-0

geno$alt_err2<-NA
geno$alt_err2_count<-0

geno$alt_err3<-NA
geno$alt_err3_count<-0


geno$alt_err1[subjectHits(error1)]<-alt_count$alt[queryHits(error1)]
geno$alt_err1_count[subjectHits(error1)]<-alt_count$N[queryHits(error1)]

geno$alt_err2[subjectHits(error2)]<-alt_count$alt[queryHits(error2)]
geno$alt_err2_count[subjectHits(error2)]<-alt_count$N[queryHits(error2)]

geno$alt_err3[subjectHits(error3)]<-alt_count$alt[queryHits(error3)]
geno$alt_err3_count[subjectHits(error3)]<-alt_count$N[queryHits(error3)]


geno$total_error<-geno$alt_err1_count+geno$alt_err2_count+geno$alt_err3_count
length(which(geno$total_error>0))
dim(geno)

pdf(file="~/plot/ASE/test.pdf", width = 10, height = 6)
ggplot(geno, aes(x=ref_ratio,y=REF_RATIO))+
  geom_point(alpha=0.3)

ggplot(geno, aes(x=ref_ratio,y=REF_RATIO))+
  geom_point(alpha=0.3)

ggplot(geno[which(geno$total_error>0),], aes(x=ref_ratio,y=REF_RATIO))+
  geom_point(alpha=0.3)
  dev.off()
  


