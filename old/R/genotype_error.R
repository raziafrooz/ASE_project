setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(ggplot2)

#sra<-readRDS("data/good_ASE.rds")

study<-"Artery_Coronary"#geno_met$study[i]
tissue<-"Artery_Coronary"   #study
geno_met<-read.csv("data/GTEx_geno_metadata.csv")
#geno<-readRDS(geno_met$allGenotypesOutput[geno_met$study==study])
ase_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rds"))
sample_id<-"GTEX-1117F-0526-SM-5EGHJ.1" #unique(ase_df$sample_id_rep)[1]
geno<- joined_filt#ase_df %>% filter(sample_id_rep==sample_id)
#ase_df<-geno %>% filter(sample_id_rep==sample_id)

#Gtex paper suggests removing vriants in HLA genes:
#his BED file contains genomic positions that we have identified as either showing bias in simulations or having a UCSC mappability score < 50. 
#Variants that fall into these positions are used for phasing, but not for generating haplotypic counts to avoid problems with mapping bias.
#https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
#Downloaded blacklist of snps from gtex (https://github.com/secastel/phaser/blob/master/phaser/README.md)
bad_snp<-read.table("~/plot/ASE/hg38_haplo_count_blacklist.chr.bed", sep="\t")
bad_snp_gr<-makeGRangesFromDataFrame(bad_snp,seqnames="V1",start.field ="V2",end.field = "V3")


library(rtracklayer)
library(GenomicRanges)
#I found a bed file including the problematic sites from ENCODE
#https://github.com/Boyle-Lab/Blacklist/tree/master
black_list<-import("~/ASE/data/hg38-blacklist.v2.bed.gz")  

#https://bismap.hoffmanlab.org
bed_graph<-import("~/k50.umap.bed")

bad_snp_gr[1:3,]
black_list[1:3]
bed_graph[1:3]



joined[1,]


joined_gr<-makeGRangesFromDataFrame(joined, seqnames.field = "chr", start.field = "start", end.field = "start")
ov<-findOverlaps(joined_gr,bad_snp_gr)
indx1<-unique(queryHits(ov))
joined_gr<-joined_gr[-unique(queryHits(ov))]

ov<-findOverlaps(joined_gr,black_list)
indx2<-unique(queryHits(ov))
joined_gr<-joined_gr[-unique(queryHits(ov))]

ov<-findOverlaps(joined_gr,bed_graph)
joined_gr$score[queryHits(ov)]<-bed_graph$score[subjectHits(ov)]
indx3<-which(joined_gr$score>=1)

joined_filt<-joined[-indx1,][-indx2,][indx3,]
joined_filt[(which(joined_filt$true_genotype!=2)),][1:4]



met_gtex<-read.csv("data/GTEx_metadata.csv")
recount3_chr_mapping <- "/dcl02/lieber/ajaffe/recount-pump/recount3.alts.chromosome_mappings.tsv"
snps_gr<-readRDS("/dcs04/hansen/data/recount_genotype/biallelic_snps/biallelic_SNP_gr.rds")
alt_path<-met_gtex$alt[which(met_gtex$sample_id_rep==sample_id & met_gtex$tissue == tissue )]
bigWig_path<-met_gtex$total[which(met_gtex$sample_id_rep==sample_id & met_gtex$tissue == tissue)]


geno_gr<-makeGRangesFromDataFrame(geno, start.field = "start", end.field = "start", keep.extra.columns = T)
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


#culate the sequencing counts aligned to non alt/ref sequence
#this should be subtracted from the total later
not_alt_count<-alt_count[-alt_key_idx,]
not_alt_count<-not_alt_count %>% group_by(chr,pos) %>% summarize(error=sum(N))
error_gr<-makeGRangesFromDataFrame(not_alt_count,seqnames="chr",start.field ="pos",end.field = "pos", keep.extra.columns = T)
colnames(not_alt_count)[2]<- "start"

geno<-right_join(not_alt_count,geno, by=c("chr","start"))

geno<-geno %>% mutate(err_per= error/(coverage.x+error))
mean(geno$err_per,na.rm=T)

colnames(geno)[7:8]<-c("coverage","ref")
geno$alt<-geno$coverage- geno$ref

norm_dat <- data.frame(q = xx$ref_ratio, cdf = yy)
ggplot(norm_dat) + geom_line(aes(x = q, y = cdf))
dev.off()

geno$alt[1:10]
geno$coverage[1:10]
geno$geno_err[1:10]
aa[1:10]

aa<-1-pbinom((geno$alt-1), size=geno$coverage, prob=mean(geno$err_per,na.rm=T))
rr<-1-pbinom((geno$ref-1), size=geno$coverage, prob=mean(geno$err_per,na.rm=T))


geno$geno_err<-rr+aa

geno<-geno %>% 
  mutate(mean=(log10(ref)+log10(alt))/2, ratio= log10(ref)-log10(alt)) %>% 
  rowwise() %>% 
  mutate(min_all=min(ref,alt))

geno_2<-geno[-which(geno$err_per>=0.05),] 

geno_final<-geno_2[-which(geno_2$geno_err>=0.001),]

pdf(file="~/plot/ASE/test/geno_error.pdf", width = 10, height = 6)

# ggplot(geno, aes(x=mean,y=ratio))+
#   geom_point(alpha=0.2)+
#   geom_point(data=geno[which(geno$true_genotype!=2),], aes(x=mean,y=ratio, color=true_genotype),alpha=0.3)
#   

ggplot(geno, aes(x=log10(ref),y=log10(alt)))+
   geom_point(alpha=0.1, color="black")+
  # geom_point(data=geno[which(geno$err_per>=0.05),], aes(x=log10(ref),y=log10(alt)),color="pink", alpha=0.3)+
  # geom_point(data=geno[which(geno$geno_err>=0.0001),], aes(x=log10(ref),y=log10(alt)),color="red", alpha=0.3)+
geom_point(data=geno[which(geno$true_genotype!=2),], aes(x=log10(ref),y=log10(alt), color=true_genotype),alpha=0.3)+
  xlim(c(0,3))+
  ylim(c(0,3))+
  labs(title="W/o any filtering. color is true geno (true_error)",
       subtitle = paste0("one gtex sample. Total snp=",nrow(geno),
                         ", total geno error=",length(which(geno$true_genotype!=2)) ))

ggplot(geno_2, aes(x=coverage,y=min_all))+
  geom_point(alpha=0.1)+
  geom_point(data=geno_2[which(geno_2$true_genotype!=2),], aes(x=coverage,y=min_all, color=true_genotype),alpha=0.1)+
  geom_point(data=geno[which(geno$err_per>=0.05),], aes(x=coverage,y=min_all),color="pink")+
  xlim(c(0,500))+
  ylim(c(0,250))+
  labs(title="error percentage in pink = error_count/total >=0.05",
       subtitle = paste0("one gtex sample. Total snp=",nrow(geno_2),
                         ", total geno error=", length(which(geno_2$true_genotype!=2)),
                         ", # snps removed here= ", sum(geno$err_per>=0.05,na.rm=T)))

ggplot(geno_2, aes(x=coverage,y=min_all))+
  geom_point(alpha=0.1)+
  geom_point(data=geno_2[which(geno_2$true_genotype!=2),], aes(x=coverage,y=min_all, color=true_genotype),alpha=0.1)+
  geom_point(data=geno_2[which(geno_2$geno_err>=0.001),], aes(x=coverage,y=min_all),color="pink")+
  xlim(c(0,500))+
  ylim(c(0,250))+
  labs(title="Tuuli method in pink(CDF) <0.001",
       subtitle = paste0("one gtex sample. Total snp=",nrow(geno_final),
                         ", total geno error=", length(which(geno_final$true_genotype!=2)),
                         ", # snps removed here= ", sum(geno_2$geno_err>=0.001,na.rm=T)))
ggplot(geno_final, aes(x=log10(ref),y=log10(alt)))+
  geom_point(alpha=0.2)+
  geom_point(data=geno_final[which(geno_final$true_genotype!=2),], aes(x=log10(ref),y=log10(alt), color=true_genotype))+
  labs(title="after both filtering",
       subtitle = paste0("one gtex sample. Total snp=",nrow(geno_final),
                         ", total geno error=", length(which(geno_final$true_genotype!=2))))+
  xlim(c(0,3))+
  ylim(c(0,3))



dev.off()



length(which(geno$geno_err>0.000001))

dim(xx)



xx<-geno[which(geno$geno_err>0.000001),] %>% filter(true_genotype!=2)



library(AnnotationHub)
ah <- AnnotationHub()
ah <- query(ah, c("v26","GENCODE","Homo sapiens","GRch38")) #v26 was used in GTExV8
TxDb<-ah[["AH75155"]]
#gene_grch38<-genes(TxDb,columns=c("TXID", "TXNAME"))
tx_38<-exons(TxDb,columns=c("TXNAME","GENEID","EXONID"))

ov<-findOverlaps(tx_38,geno_gr)

geno_gr$exon<-"no"
geno_gr$exon[subjectHits(ov)]<-"yes"
table(geno_gr$exon,geno_gr$true_genotype)


xx<-geno[which(geno$err_per>=0.05),]
xx<-geno_2[which(geno_2$geno_err>=0.001),]
geno[1,]
length(which(xx$true_genotype!=2))
dim(xx)
length(which(xx$true_genotype!=2))/dim(xx)[1]
