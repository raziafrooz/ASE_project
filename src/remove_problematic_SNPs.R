library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library("org.Hs.eg.db")


#Load in the problematic sites from publicly available dataset:
#---------------------------------------
#Downloaded blacklist of snps from gtex (https://github.com/secastel/phaser/blob/master/phaser/README.md)
HLA_snp<-read.table("/users/arazi/hansen_lab/dwl_files/ASE_filter/hg38_haplo_count_blacklist.chr.bed", sep="\t")
HLA_snp_gr<-makeGRangesFromDataFrame(HLA_snp,seqnames="V1",start.field ="V2",end.field = "V3")

#I found a bed file including the problematic sites from ENCODE
#https://github.com/Boyle-Lab/Blacklist/tree/master
black_list<-import("/users/arazi/hansen_lab/dwl_files/ASE_filter/hg38-blacklist.v2.bed.gz")

#https://bismap.hoffmanlab.org
mappability<-import("/users/arazi/hansen_lab/dwl_files/ASE_filter/k50.Unique.Mappability.bb")


#----------------------------------------
#GTEX has some sites marked as bad. load these in an remove it from our data:
gtex_sim<-fread("/dcs07/hansen/data/recount_ASE/data/gtex_simulation.csv.gz")
gtex_sim_gr<-makeGRangesFromDataFrame(gtex_sim,seqnames="chr",start.field ="start",end.field = "start")
#===========================================================
ah <- AnnotationHub()
ah <- query(ah, c("v26","GENCODE","Homo sapiens","GRch38")) #v26 was used in GTExV8
TxDb<-ah[["AH75155"]]
tx_38<-exons(TxDb,columns=c("TXNAME","GENEID","EXONID","CDSID"))
#===========================================================




remove_problematic_SNPs<-function(ase_df,bigWig_path){
  
  #bigWig_path<-gtex_predict$total[which(gtex_predict$sample_id_rep==sample_id)][1]
  ase_df_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")
  
  
  
  bigwig <-tryCatch(
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
  overlap_loci <- findOverlaps(ase_df_gr, bigwig)
  ase_df$bigwig_count <- bigwig$score[subjectHits(overlap_loci)]
  
  #Make Granges:
  ase_filt_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")
  
  #Remove blacklist snps from our granges:
  ov<-findOverlaps(ase_filt_gr,HLA_snp_gr)
  ase_df<-ase_df[-unique(queryHits(ov)),]
  ase_filt_gr<-ase_filt_gr[-unique(queryHits(ov)),]
  n1<-length(unique(queryHits(ov)))
  
  ov<-findOverlaps(ase_filt_gr,black_list)
  ase_df<-ase_df[-unique(queryHits(ov)),]
  ase_filt_gr<-ase_filt_gr[-unique(queryHits(ov)),]
  n1<-length(unique(queryHits(ov)))+n1
  
  ov<-findOverlaps(ase_filt_gr,mappability)
  ase_df<-ase_df[unique(queryHits(ov)),]
  ase_df_gr<-ase_filt_gr[unique(queryHits(ov)),]
  n1<-length(ase_filt_gr[-unique(queryHits(ov)),])+n1
  
  print(paste0("The number of SNPs removed based on publicly available data is = ",n1))
  #ge
  #---------------------
  
  ov<-findOverlaps(ase_df_gr,gtex_sim_gr)
  ase_df$LOW_MAPABILITY<-NA
  ase_df$LOW_MAPABILITY[queryHits(ov)]<-gtex_sim$LOW_MAPABILITY[subjectHits(ov)]
  ase_df$MAPPING_BIAS_SIM<-NA
  ase_df$MAPPING_BIAS_SIM[queryHits(ov)]<-gtex_sim$MAPPING_BIAS_SIM[subjectHits(ov)]
  
  n2<-nrow(ase_df %>% filter(LOW_MAPABILITY>0,MAPPING_BIAS_SIM>0))
  
  ase_df<-ase_df %>% filter(LOW_MAPABILITY<1,MAPPING_BIAS_SIM<1)
  print(paste0("The number of SNPs removed based on GTEx simulation is = ",n2))
  
  #------
  #remove HLA GENES
  
  
  ase_df_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")
  ov<-findOverlaps(tx_38,ase_df_gr)
  
  ase_df$GENE_ID<-NA
  ase_df$GENE_ID[subjectHits(ov)]<-unlist(tx_38$GENEID[queryHits(ov)])
  ase_df$GENE_ID<-gsub("\\.\\d+", "", ase_df$GENE_ID)
  
  gene_sym<- mapIds(org.Hs.eg.db,
                    keys=ase_df$GENE_ID, #Column containing Ensembl gene ids
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
  id<-match(ase_df$GENE_ID, names(gene_sym))
  gene_sym<-gene_sym[!is.na(names(gene_sym))]
  gene_sym<-data.frame(symbol=unlist(gene_sym), GENE_ID=names(unlist(gene_sym)))
  
  ase_df$symbol<-gene_sym$symbol[match(ase_df$GENE_ID,gene_sym$GENE_ID)]
  
  n3<-nrow(ase_df %>%filter(str_detect(symbol,"HLA")))
  
  ase_df<-ase_df %>%
    filter(!str_detect(symbol,"HLA"))
  
  print(paste0("The number of remaining HLA genes removed is = ",n3))
  
    return(ase_df)
}
  
  
  
  
  
  
  
  
  


