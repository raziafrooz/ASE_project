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



gtex_w_ov<-qc_df %>% filter(SMGEBTCHT=="TruSeq.v1",overlap> 40) %>% 
  left_join(gtex_uni_norm,by=c("SAMPLE_ID")) %>% 
  filter(uni_norm_mean<0.08)
row_id=2
recount_sig_snps<-c()
for(row_id in 1:nrow(gtex_w_ov)){
  sample_id<-gtex_w_ov$SAMPLE_ID[row_id]
  study<-gtex_w_ov$study[row_id]
  overlap<-gtex_w_ov$overlap[row_id]
  indv<-gtex_w_ov$SUBJID[row_id]
  print(row_id)
  
  ase_df<-fread(gtex_metadata$genotypedSamples[which(gtex_metadata$sample_id_rep==sample_id)][1]) %>% 
    filter(pred_genotype==2, coverage>=8) %>% 
    mutate(ref_ratio=ref_count/coverage,
           alt_count=coverage-ref_count,
           ratio=log2(ref_count/alt_count),
           mean=(log2(ref_count)+log2(alt_count))/2)
  
  #Fix the counts by adjusting the MA plot to mean around 0
  ratio_adj<-median(ase_df$ratio)
  if(ratio_adj!=0){
    
    adj<-ase_df$ratio-ratio_adj
    ase_df$adj_alt<-round(ase_df$coverage/((2^adj)+1),0)
    
    stopifnot(sum(is.na(ase_df$adj_alt))==0)
    
    ase_df$adj_ref<-ase_df$coverage-ase_df$adj_alt
    
  }
  
  
  
  
  
  
  
  
  bigWig_path<-xx$total[which(xx$sample_id==sample_id)]
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
    dplyr::select(chr, start,ref_count,alt_count, coverage,bigwig_count,ref_ratio,adj_alt,adj_ref )
  
  
  ase_df$err_per <- (ase_df$bigwig_count - ase_df$coverage)/ase_df$bigwig_count
  
  
  aa<-1-pbinom((ase_df$alt_count-1), size=ase_df$coverage, prob=mean(ase_df$err_per,na.rm=T))
  rr<-1-pbinom((ase_df$ref_count-1), size=ase_df$coverage, prob=mean(ase_df$err_per,na.rm=T))
  
  
  ase_df$geno_err<-rr+aa
  
  
  ase_df<-ase_df[-which(ase_df$err_per>=0.05),] 
  ase_df<-ase_df[-which(ase_df$geno_err>=0.001),]
  
  
  
  
  ase_df_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")
  
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
  
  
  
  
  is.median<-median(ase_df$adj_ref/ase_df$coverage)
  ase_df$p_val = apply(ase_df[,c("adj_ref","adj_alt")], 1, function(x) {
    binom.test(round(x[1],1),round((x[1]+x[2]),1),p=is.median)$p.value})
  
  ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")
  
  wasp<-tissues_names$file_name[tissues_names$indv==indv]
  if(length(wasp)>0){
    gtex_ase<- fread(wasp)
    colnames(gtex_ase)[1:2]<- c("chr", "start")
    gtex_ase_1<-gtex_ase %>% 
      filter(SAMPLE_ID ==sample_id,LOW_MAPABILITY<1, MAPPING_BIAS_SIM<1, GENOTYPE_WARNING<1) %>%
      dplyr::select(c(chr, start,REF_COUNT,ALT_COUNT,TOTAL_COUNT,REF_RATIO,BINOM_P_ADJUSTED,BINOM_P))
    
    #         joined<-ase_df  %>% 
    #           inner_join(gtex_ase_1,by = c("chr", "start")) 
  }
  
  dd<-ase_df  %>% full_join(gtex_ase_1,by = c("chr", "start"))
  dd_union<-inner_join(gtex_ase_1,ase_df,by=c("chr","start"))
  
  n_union<-nrow(dd_union)
  n_recount<-sum(!is.na(dd$ref_count))
  n_gtex<-sum(!is.na(dd$REF_COUNT))
  #before union:
  sig_recount.05=sum(dd$q_val<0.05, na.rm=T)
  sig_gtex= sum(dd$BINOM_P_ADJUSTED<0.05,na.rm=T)
  gene_rec<-sum(!is.na((unique(dd$GENE_ID[dd$q_val<0.05]))))
  gene_wasp<-sum(!is.na(unique(dd$GENE_ID[dd$BINOM_P_ADJUSTED<0.05])))
  gene_both<-sum(!is.na(unique(dd$GENE_ID[dd$BINOM_P_ADJUSTED<0.05 & dd$q_val<0.05])))
  
  
  sig_recount.05_uni=sum(dd_union$q_val<0.05, na.rm=T)
  sig_gtex_uni= sum(dd_union$BINOM_P_ADJUSTED<0.05,na.rm=T)
  sig_both_uni= sum(dd_union$BINOM_P_ADJUSTED<0.05 & dd_union$q_val<0.05,na.rm=T)
  uni_gene_rec<-sum(!is.na((unique(dd_union$GENE_ID[dd_union$q_val<0.05]))))
  uni_gene_wasp<-sum(!is.na(unique(dd_union$GENE_ID[dd_union$BINOM_P_ADJUSTED<0.05])))
  
  df_x<-data.frame(SAMPLE_ID=sample_id,
                   overlap,
                   total_n_intersect=n_union,
                   total_n_recount=n_recount,
                   total_n_gtex=n_gtex,
                   snp_sig_onlyRecount=sig_recount.05,
                   snp_sig_onlyGtex=sig_gtex,
                   snp_sig_Recount=sig_recount.05_uni,
                   snp_sig_gtex=sig_gtex_uni,
                   snp_sig_both=sig_both_uni,
                   gene_sig_onlyRecount=gene_rec,
                   gene_sig_onlyGtex=gene_wasp,
                   gene_sig_Recount=uni_gene_rec,
                   gene_sig_Gtex=uni_gene_wasp, gene_both,study)
  
  recount_sig_snps<-rbind(recount_sig_snps,df_x)
}
}
}

plot_df2<-ase_df  %>% full_join(gtex_ase_1,by = c("chr", "start"))
plot_df<-ase_df  %>% full_join(gtex_ase_1,by = c("chr", "start")) %>% 
  filter(BINOM_P_ADJUSTED>=0.05 & q_val<0.05)

plot_df3<-ase_df  %>% full_join(gtex_ase_1,by = c("chr", "start")) %>% 
  filter(BINOM_P_ADJUSTED<0.05 & q_val>=0.05)

pdf(file="~/plot/ASE/high_ov2.pdf", width = 10, height = 6)
ggplot(plot_df, aes(x=adj_ref/coverage , y=REF_COUNT/TOTAL_COUNT))+
  geom_point(alpha=0.3)+
  geom_point(data=plot_df, aes(x=ref_count/coverage , y=REF_COUNT/TOTAL_COUNT),alpha=0.3,color="brown")+
  geom_abline(slope=1, color="red")+
  labs(title=paste0(sample_id," SNPs only sig in recount"), 
       subtitle=paste0("brown is after adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")

ggplot(plot_df3, aes(x=adj_ref/coverage , y=REF_COUNT/TOTAL_COUNT))+
  geom_point(alpha=0.3)+
  geom_point(data=plot_df3, aes(x=ref_count/coverage , y=REF_COUNT/TOTAL_COUNT),alpha=0.3,color="brown")+
  geom_abline(slope=1, color="red")+
  labs(title=paste0(sample_id,": SNPs only sig in GTEx"), 
       subtitle=paste0("brown is after adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")

ggplot(plot_df, aes(BINOM_P))+
  geom_histogram(alpha=0.3,fill="darkgreen")+
  geom_histogram(data=plot_df, aes(p_val),alpha=0.3)+
  labs(title=paste0(sample_id,": all SNPs (green is gtex)"),
       subtitle=paste0("overlap is: ",overlap),
       x="p_value")

ggplot(plot_df, aes(REF_COUNT/TOTAL_COUNT))+
  geom_histogram(alpha=0.3,fill="darkgreen")+
  geom_histogram(data=plot_df, aes(adj_ref/coverage),alpha=0.3)+
  labs(title=paste0(sample_id,": all SNPs (green is gtex)"),
       subtitle=paste0("overlap is: ",overlap),
       x="ref_ratio")

ggplot(plot_df3, aes(BINOM_P))+
  geom_histogram(alpha=0.3,fill="darkgreen")+
  geom_histogram(data=plot_df3, aes(p_val),alpha=0.3)+
  labs(title=paste0(sample_id,": SNPs only sig in gtex (green is gtex)"),
       subtitle=paste0("overlap is: ",overlap),
       x="p_value")

ggplot(plot_df3, aes(REF_COUNT/TOTAL_COUNT))+
  geom_histogram(alpha=0.3,fill="darkgreen")+
  geom_histogram(data=plot_df3, aes(adj_ref/coverage),alpha=0.3)+
  labs(title=paste0(sample_id,": SNPs only sig in gtex  (green is gtex)"),
       subtitle=paste0("overlap is: ",overlap),
       x="ref_ratio")


ggplot(plot_df2, aes(BINOM_P))+
  geom_histogram(alpha=0.3,fill="darkgreen")+
  geom_histogram(data=plot_df2, aes(p_val),alpha=0.3)+
  labs(title=paste0(sample_id,": SNPs only sig in recount  (green is gtex)"),
       subtitle=paste0("overlap is: ",overlap),
       x="P-value")

ggplot(plot_df2, aes(REF_COUNT/TOTAL_COUNT))+
  geom_histogram(alpha=0.3,fill="darkgreen")+
  geom_histogram(data=plot_df2, aes(adj_ref/coverage),alpha=0.3)+
  labs(title=paste0(sample_id,": SNPs only sig in recount  (green is gtex)"),
       subtitle=paste0("overlap is: ",overlap),
       x="ref_ratio")

dev.off()

plot_df2[1,]
plot_df4<-plot_df2 %>% rowwise() %>%  mutate(min_allele=min(ref_count, alt_count))



pdf(file="~/plot/ASE/test3.pdf", width = 10, height = 6)
ggplot(plot_df4, aes(y=min_allele , x=coverage))+
  geom_point(alpha=0.3)+xlim(c(0,800))+ylim(c(0,400))+
  geom_point(data=plot_df4 %>% filter(BINOM_P_ADJUSTED>=0.05 & q_val<0.05),aes(y=min_allele , x=coverage),color="red")
dev.off()


#---------------------------------
plot_df<-plot_df2 %>%  mutate(ratio_total=log2(coverage/TOTAL_COUNT),
                              mean_total=(log2(coverage)+log2(TOTAL_COUNT))/2,
                              ratio_ref=log2(ref_count/REF_COUNT),
                              mean_ref=(log2(ref_count)+log2(REF_COUNT))/2,
                              ratio_alt=log2(alt_count/ALT_COUNT),
                              mean_alt=(log2(alt_count)+log2(ALT_COUNT))/2,
                              ratio_ref_adj=log2(adj_ref/REF_COUNT),
                              mean_ref_adj=(log2(adj_ref)+log2(REF_COUNT))/2,
                              ratio_alt_adj=log2(adj_alt/ALT_COUNT),
                              mean_alt_adj=(log2(adj_alt)+log2(ALT_COUNT))/2,
                              ratio=log2(ref_count/alt_count),
                              mean=(log2(ref_count)+log2(alt_count))/2,
                              ratio2=log2(adj_ref/adj_alt),
                              mean2=(log2(adj_ref)+log2(adj_alt))/2)
pdf(file="~/plot/ASE/high_ov_test2.pdf", width = 10, height = 6)
ggplot(plot_df, aes(x=mean , y=ratio))+
  geom_point(alpha=0.3)+
  labs(title=paste0(sample_id," SNPs only sig in recount"), 
       subtitle=paste0("before adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")+geom_smooth()

ggplot(plot_df, aes(x=mean2 , y=ratio2))+
  geom_point(alpha=0.3)+
  labs(title=paste0(sample_id," SNPs only sig in recount"), 
       subtitle=paste0("before adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")+geom_smooth()

ggplot(plot_df, aes(x=mean_total , y=ratio_total))+
  geom_point(alpha=0.3)+
  labs(title=paste0(sample_id," SNPs only sig in recount"), 
       subtitle=paste0("before adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")+geom_smooth()

ggplot(plot_df, aes(x=mean_ref , y=ratio_ref))+
  geom_point(alpha=0.3)+
  labs(title=paste0(sample_id," SNPs only sig in recount"), 
       subtitle=paste0("before adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")+geom_smooth()
ggplot(plot_df, aes(x=mean_ref_adj , y=ratio_ref_adj))+
  geom_point(alpha=0.3)+
  labs(title=paste0(sample_id," SNPs only sig in recount"), 
       subtitle=paste0("before adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")+geom_smooth()
ggplot(plot_df, aes(x=mean_alt , y=ratio_alt))+
  geom_point(alpha=0.3)+
  labs(title=paste0(sample_id," SNPs only sig in recount"), 
       subtitle=paste0("before adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")+geom_smooth()
ggplot(plot_df, aes(x=mean_alt_adj , y=ratio_alt_adj))+
  geom_point(alpha=0.3)+
  labs(title=paste0(sample_id," SNPs only sig in recount"), 
       subtitle=paste0("before adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")+geom_smooth()

dev.off()

ggplot(plot_df, aes(x=adj_ref/coverage , y=REF_COUNT/TOTAL_COUNT))+
  geom_point(alpha=0.3)+
  geom_abline(slope=1, color="red")+
  labs(title=paste0(sample_id," SNPs only sig in recount"), 
       subtitle=paste0("before adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")


ggplot(plot_df, aes(x=ref_count/coverage , y=REF_COUNT/TOTAL_COUNT))+
  geom_point(alpha=0.3,color="brown")+
  geom_abline(slope=1, color="red")+
  labs(title=paste0(sample_id," SNPs only sig in recount"), 
       subtitle=paste0("after adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")

dev.off()
ggplot(plot_df3, aes(x=adj_ref/coverage , y=REF_COUNT/TOTAL_COUNT))+
  geom_point(alpha=0.3)+
  geom_point(data=plot_df3, aes(x=ref_count/coverage , y=REF_COUNT/TOTAL_COUNT),alpha=0.3,color="brown")+
  geom_abline(slope=1, color="red")+
  labs(title=paste0(sample_id,": SNPs only sig in GTEx"), 
       subtitle=paste0("brown is after adj. overlap is: ",overlap),
       x="recount_refRatio",
       y="Gtex_refRatio")

