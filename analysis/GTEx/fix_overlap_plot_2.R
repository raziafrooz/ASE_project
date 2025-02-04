#-----------------
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library("org.Hs.eg.db")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
#Downloaded blacklist of snps from gtex (https://github.com/secastel/phaser/blob/master/phaser/README.md)
HLA_snp<-read.table("/users/arazi/hansen_lab/dwl_files/ASE_filter/hg38_haplo_count_blacklist.chr.bed", sep="\t")
HLA_snp_gr<-makeGRangesFromDataFrame(HLA_snp,seqnames="V1",start.field ="V2",end.field = "V3")

#I found a bed file including the problematic sites from ENCODE
#https://github.com/Boyle-Lab/Blacklist/tree/master
black_list<-import("/users/arazi/hansen_lab/dwl_files/ASE_filter/hg38-blacklist.v2.bed.gz")

#https://bismap.hoffmanlab.org
mappability<-import("/users/arazi/hansen_lab/dwl_files/ASE_filter/k50.Unique.Mappability.bb")


#===========================================================
ah <- AnnotationHub()
ah <- query(ah, c("v26","GENCODE","Homo sapiens","GRch38")) #v26 was used in GTExV8
TxDb<-ah[["AH75155"]]
tx_38<-exons(TxDb,columns=c("TXNAME","GENEID","EXONID","CDSID"))

xx<-read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/metadata/all_gtex_metadata.csv")
gtex_sim<-fread("/dcs07/hansen/data/recount_ASE/data/gtex_simulation.csv.gz")
gtex_sim_gr<-makeGRangesFromDataFrame(gtex_sim,seqnames="chr",start.field ="start",end.field = "start")


gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)

qc_df<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/gtex_qc_metadata.csv.gz")
qc_df<-qc_df %>% mutate(overlap= (star.average_input_read_length)-bc_frag.mode_length)
gtex_uni_norm<-fread("/dcs07/hansen/data/recount_ASE/data/gtex_uniNorm.csv")
colnames(gtex_uni_norm)[2]<-"SAMPLE_ID"
qc_df$SAMPLE_ID<-str_sub(qc_df$external_id, end= -3)

#============================================================
#Read in wasp indv:
path_id<-"/dcl01/hansen/data/arazi/ASE/dbGap/phe000039.v1.GTEx_v8_ASE_WASP.expression-matrixfmt-ase.c1/GTEx_Analysis_v8_ASE_WASP_counts_by_subject/"
tissues_names<-as.data.frame(list.files(path =path_id))
colnames(tissues_names)<-"file_name"
tissues_names$indv<-gsub("\\..*","",tissues_names$file_name)
tissues_names$file_name<-paste0(path_id, tissues_names$file_name)
#============================================================

gtex_w_ov<-qc_df %>% filter(SMGEBTCHT=="TruSeq.v1") %>% 
  left_join(gtex_uni_norm,by=c("SAMPLE_ID")) %>% 
  filter(uni_norm_mean<0.08)
gtex_w_ov[which(gtex_w_ov$overlap> 50),]
adj_stat_all[1:3,]
#--------------------------
#adj_stat_all<-c()
row_id=which(gtex_w_ov$SAMPLE_ID=="GTEX-NPJ8-1526-SM-2D7VU")
for(row_id in 1:nrow(gtex_w_ov)){
  sample_id<-gtex_w_ov$SAMPLE_ID[row_id]
  study<-gtex_w_ov$study[row_id]
  overlap<-gtex_w_ov$overlap[row_id]
  indv<-gtex_w_ov$SUBJID[row_id]
  print(row_id)
  saveTempFile<-paste0("~/hansen_lab/ASE/test/GTEx/",sample_id,".csv.gz")
  if(!file.exists(saveTempFile)){
  ase_df<-fread(gtex_metadata$genotypedSamples[which(gtex_metadata$sample_id_rep==sample_id)][1]) %>% 
    filter(pred_genotype==2, coverage>=8) %>% 
    mutate(ref_ratio=ref_count/coverage,
           alt_count=coverage-ref_count,
           ratio=log2(ref_count/alt_count),
           mean=(log2(ref_count)+log2(alt_count))/2)
  
  wasp<-tissues_names$file_name[tissues_names$indv==indv]
  if(length(wasp)>0){
  gtex_ase<- fread(wasp)
  colnames(gtex_ase)[1:2]<- c("chr", "start")
  gtex_ase_1<-gtex_ase %>% 
    filter(SAMPLE_ID ==sample_id,LOW_MAPABILITY<1, MAPPING_BIAS_SIM<1, GENOTYPE_WARNING<1) #%>%
    #dplyr::select(c(chr, start,REF_COUNT,ALT_COUNT,TOTAL_COUNT,REF_RATIO,BINOM_P_ADJUSTED,BINOM_P))
  
  if(length(gtex_ase_1)>0){
  
  
  bigWig_path<-xx$total[which(xx$sample_id==sample_id)][1]
  #bigWig_path<-xx_one$total[which(xx_one$sample_id_rep==sam_id)]
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
  
  
  
  ase_df<-ase_df %>%
    dplyr::select(chr, start,ref_count,alt_count, coverage,bigwig_count,ref_ratio,ratio,mean )


  ase_df$err_per <- (ase_df$bigwig_count - ase_df$coverage)/ase_df$bigwig_count
  error_prob=mean(ase_df$err_per,na.rm=T)
  
  
  aa<-1-sapply(1:nrow(ase_df), function(zz)  { pbinom((ase_df$alt_count[zz]-1), size=ase_df$coverage[zz], prob=error_prob) })
  rr<-1-sapply(1:nrow(ase_df), function(zz)  { pbinom((ase_df$ref_count[zz]-1), size=ase_df$coverage[zz], prob=error_prob) })
  
  
  
  ase_df$geno_err<-rr+aa
  
  
  ase_df<-ase_df[-which(ase_df$err_per>=0.05),] 
  ase_df<-ase_df[-which(ase_df$geno_err>=0.001),]
  
  #----
  
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
  ase_df_gr<-ase_filt_gr[unique(queryHits(ov)),]
  #ge
  #---------------------
  
  
  #ase_df_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")
  
  ov<-findOverlaps(ase_df_gr,gtex_sim_gr)
  ase_df$LOW_MAPABILITY<-NA
  ase_df$LOW_MAPABILITY[queryHits(ov)]<-gtex_sim$LOW_MAPABILITY[subjectHits(ov)]
  ase_df$MAPPING_BIAS_SIM<-NA
  ase_df$MAPPING_BIAS_SIM[queryHits(ov)]<-gtex_sim$MAPPING_BIAS_SIM[subjectHits(ov)]
  
  ase_df<-ase_df %>% filter(LOW_MAPABILITY<1,MAPPING_BIAS_SIM<1)
  
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
  
  ase_df<-ase_df %>%
    filter(!str_detect(symbol,"HLA"))
  
  #------
  fwrite(ase_df,saveTempFile)
  }
  else{
    ase_df<-fread(saveTempFile)
  }
  #Fix the counts by adjusting the MA plot to mean around 0
  ratio_adj<-median(ase_df$ratio)
  if(ratio_adj!=0){
    
  
  adj<-ase_df$ratio-ratio_adj
  ase_df$adj_alt<-round((ase_df$coverage/((2^adj)+1)),0)
  
  stopifnot(sum(is.na(ase_df$adj_alt))==0)

  ase_df$adj_ref<-ase_df$coverage-ase_df$adj_alt
  
  is.median<-median(ase_df$adj_ref/ase_df$coverage)
  ase_df$p_val = apply(ase_df[,c("adj_ref","adj_alt")], 1, function(x) {
    binom.test(round(x[1],1),round((x[1]+x[2]),1),p=is.median)$p.value})
  
  ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")
  

    #         joined<-ase_df  %>% 
    #           inner_join(gtex_ase_1,by = c("chr", "start")) 
  
  
  joined<-ase_df %>% 
    dplyr::select(c(chr, start,ref_count,alt_count,adj_alt, adj_ref,coverage,ref_count,alt_count,p_val,q_val)) %>% 
    inner_join(gtex_ase_1,by = c("chr", "start")) 
  
  
  adj_stat<-data.frame(overlap,
                       sample_id,
                       ref_before_adj=median((joined$ref_count -joined$REF_COUNT)/joined$ref_count),
                       ref_after_adj=median((joined$adj_ref -joined$REF_COUNT)/joined$adj_ref),
                       alt_before_adj=median((joined$alt_count -joined$ALT_COUNT)/joined$alt_count),
                       alt_after_adj=median((joined$adj_alt -joined$ALT_COUNT)/joined$adj_alt),
                       total_diff=median((joined$coverage -joined$TOTAL_COUNT)/joined$coverage))
  
  
  adj_stat_all<-rbind(adj_stat_all,adj_stat)
  
  } 
}
  }
}
if(!file.exists("~/test/gtex_fixed.csv")){fwrite(adj_stat_all, file="~/test/gtex_fixed.csv")}
adj_stat_all<-fread("~/test/gtex_fixed.csv")

#--------------------
#Make some plots:
pp_df<-adj_stat_all %>%
  mutate(overlap_sign=case_when(overlap >= 0 ~ "Positive",
                                overlap < 0 ~ "Negative"))  %>% 
  group_by(overlap_sign) %>% 
  mutate(ov_group=cut_number(overlap,3))
pp_df_long<- pp_df %>%
  dplyr::select(overlap_sign,overlap,sample_id,ov_group,alt_before_adj,alt_after_adj,total_diff) %>% 
  pivot_longer(!c(overlap_sign,overlap,sample_id,ov_group), names_to="adjustment", values_to="difference")

pdf(file="~/plot/ASE/adjustVSwasp_test3.pdf", width = 10, height = 6)

pp_df_long<- pp_df %>% dplyr::select(overlap_sign,overlap,sample_id,ov_group,ref_before_adj,alt_before_adj,total_diff) %>% pivot_longer(!c(overlap_sign,overlap,sample_id,ov_group), names_to="adjustment", values_to="difference")

ggplot(data=pp_df_long ,aes(y=difference, x=adjustment, color=ov_group))+
  geom_boxplot(alpha=0.4,outliers=F)+
  ylim(c(0,0.15))+
  facet_wrap(vars(overlap_sign))+
  labs(title=paste0(nrow(pp_df)," samples before adjustments"))

pp_df_long<- pp_df %>% dplyr::select(overlap_sign,overlap,sample_id,ov_group,ref_after_adj,alt_after_adj,total_diff) %>% pivot_longer(!c(overlap_sign,overlap,sample_id,ov_group), names_to="adjustment", values_to="difference")

ggplot(data=pp_df_long ,aes(y=difference, x=adjustment, color=ov_group))+
  geom_boxplot(alpha=0.4,outliers=F)+
  ylim(c(0,0.15))+
  facet_wrap(vars(overlap_sign))+
  labs(title=paste0(nrow(pp_df)," samples after adjustments"))

ggplot(data=pp_df %>% filter(alt_after_adj>=0),aes(y=alt_after_adj, x=ref_after_adj, color=ov_group))+
  geom_point(alpha=0.4)+
  labs(title="After correcting overlap issue in 4,367 GTEx samples")+
  geom_abline(slope=1, color="red")+
  facet_wrap(vars(overlap_sign))+
  labs(title=paste0(nrow(pp_df)," samples after adjustments"))

ggplot(data=pp_df %>% filter(alt_after_adj>=0),aes(y=alt_before_adj, x=ref_before_adj, color=ov_group))+
  geom_point(alpha=0.4)+
  labs(title="before correcting overlap issue in 4,367 GTEx samples")+
  geom_abline(slope=1, color="red")+
  facet_wrap(vars(overlap_sign))+
  labs(title=paste0(nrow(pp_df)," samples before adjustments"))
dev.off()





pdf(file="~/plot/ASE/atest3.pdf", width = 10, height = 6)

pp_df<-adj_stat_all %>%
  mutate(overlap_sign=case_when(overlap >= 0 ~ "Positive",
                                overlap < 0 ~ "Negative"))  %>% 
  group_by(overlap_sign) %>% 
  mutate(ov_group=cut_number(overlap,3))

ggplot(data=pp_df ,aes(y=total_diff, x=ov_group, color=overlap_sign))+
  geom_boxplot(alpha=0.4,outliers=F)

pp_df_long<- pp_df %>% dplyr::select(overlap_sign,overlap,sample_id,ov_group,alt_before_adj,alt_after_adj) %>% pivot_longer(!c(overlap_sign,overlap,sample_id,ov_group), names_to="adjustment", values_to="difference")

ggplot(data=pp_df_long ,aes(y=difference, x=ov_group, color=adjustment))+
  geom_boxplot(alpha=0.4,outliers=F)
# ggplot(data=pp_df ,aes(y=alt_after_adj, x=ov_group, color=overlap_sign))+
#   geom_boxplot(alpha=0.4,outliers=F)
dev.off()

# +
#   ylim(c(0,0.15))+
#   facet_wrap(vars(overlap_sign))+
#   labs(title=paste0(nrow(pp_df)," samples before adjustments"))
# 
# 
# 
# 
