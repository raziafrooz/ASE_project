#------------------------------------------
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
ah <- AnnotationHub()
ah <- query(ah, c("v26","GENCODE","Homo sapiens","GRch38")) #v26 was used in GTExV8
TxDb<-ah[["AH75155"]]

#Downloaded blacklist of snps from gtex (https://github.com/secastel/phaser/blob/master/phaser/README.md)
HLA_snp<-read.table("/users/arazi/hansen_lab/dwl_files/ASE_filter/hg38_haplo_count_blacklist.chr.bed", sep="\t")
HLA_snp_gr<-makeGRangesFromDataFrame(HLA_snp,seqnames="V1",start.field ="V2",end.field = "V3")

#I found a bed file including the problematic sites from ENCODE
#https://github.com/Boyle-Lab/Blacklist/tree/master
black_list<-import("/users/arazi/hansen_lab/dwl_files/ASE_filter/hg38-blacklist.v2.bed.gz")

#https://bismap.hoffmanlab.org
mappability<-import("/users/arazi/hansen_lab/dwl_files/ASE_filter/k50.Unique.Mappability.bb")
#-----------------
all_metadata<-fread("/dcs07/hansen/data/recount_ASE/data/tcga_recount_metadata.csv")
tcga<-read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/tcga.normal.metadata.csv")
raw_df<-read_csv("~/recount_genotype/redo_TCGA_snakemake/tcga.normal.metadata.csv")

sig_all<-c()
median_ref<-c()


for(i in 1:nrow(tcga)){
  print(i)
  ase_df<-fread(tcga$genotypedSamples[i])
  
  ase_df<-ase_df %>% 
    filter(coverage>=8, pred_genotype==2) %>% 
    mutate(alt_count=coverage-ref_count,
           ref_ratio=ref_count/coverage,
           log2aFC= log2((alt_count+1)/(ref_count+1)))
  
  median(ase_df$ref_ratio)
  
  
  
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
  #ge
  #-------------------------------
  ase_df_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")
  
  tx_38<-exons(TxDb,columns=c("TXNAME","GENEID","EXONID","CDSID"))
  
  ov<-findOverlaps(tx_38,ase_df_gr)
  
  ase_df$GENEID<-NA
  ase_df$GENEID[subjectHits(ov)]<-unlist(tx_38$GENEID[queryHits(ov)])
  
  
  
  ase_df<-ase_df %>% 
    dplyr::select(-c(AF,M,S,predicted_accuracy,pred_genotype,chr,start)) %>% 
    group_by(GENEID) %>% 
    summarise(coverage=round(mean(coverage),0),
              ref_count=round(mean(ref_count),0),
              alt_count=round(mean(alt_count),0),
              ref_ratio=round(mean(ref_ratio),3),
              log2aFC=round(mean(log2aFC),2)) %>% 
    ungroup()
  #-------------------------------
  
  
  is.median<-round(median(ase_df$ref_ratio),3)
  
  print(paste0("median is ",is.median))
  # Get the binomial p-value
  ase_df$p_val = apply(ase_df[,c("ref_count","alt_count")], 1, function(x) {
    binom.test(round(x[1],1),round((x[1]+x[2]),1),p=is.median)$p.value})
  
  # perform multiple testing correction with FDR
  ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")
  sig<-ase_df %>% filter(q_val<0.05) %>% 
    dplyr::select(GENEID,log2aFC)
  colnames(sig)[2]<-paste0("log2aFC", "_",i)
  
  if(i==1){
    sig_all<-sig
  }else{
    sig_all<-full_join(sig_all,sig,by='GENEID')
  }
  median_ref<-rbind(median_ref,is.median)
}

median_df<-as.data.frame(median_ref)
median_df$sample_id<- tcga$sample_id
colnames(median_df)[1]<-"ref_ratio"

# fwrite(sig_all,"~/test/tcga.csv")
# fwrite(median_df,"~/test/tcga_ref_ratio.csv")
# sig_all<-fread("~/test/tcga.csv")
# median_ref<-fread("~/test/tcga_ref_ratio.csv")



sig_all[1:10,1:6]

median_df[1,]
all_metadata[1,1:10]
str_detect(). colnames(all_metadata)
colnames(all_metadata)[str_detect(colnames(all_metadata), "mode_length" )]

median_df$bc_frag.mode_length<-all_metadata$bc_frag.mode_length[match(median_df$sample_id,all_metadata$tcga_barcode)]
median_df$star.average_input_read_length<-as.numeric(all_metadata$star.average_input_read_length[match(median_df$sample_id,all_metadata$tcga_barcode)])


median_df[,c("ref_ratio","overlap")][1:10,]

median_df<-median_df %>% mutate(overlap= star.average_input_read_length-bc_frag.mode_length)

#Single end should  ot have bc_frag.mode_length, but some have >0 values--> fix them
ase_df$overlap[which(ase_df$library_layout=="single" & ase_df$bc_frag.mode_length==0)]<-ase_df$star.average_input_read_length[which(ase_df$library_layout=="single" & ase_df$bc_frag.mode_length==0)]

sig_all[which(id<300),][1:4,1:5]
apply(sig_all,1,sum(is.na()))

sig_all %>% rowwise() %>% summarize(is.na())
sig_all[1:4,1:13]


id<-rowSums(is.na(sig_all))

quantile(id)





