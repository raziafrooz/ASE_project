library(tidyverse)
library(data.table)
library(VGAM)

#Define the functions
#-------------------------------------------------
compare_binom <- function(coverage,ref_count){
  sizes = coverage
  random = sapply(sizes, function(zz) rbinom(size = zz, n = 1, p = 0.5))
  suppressWarnings({ ks.test(random/sizes, ref_count/sizes)$statistic})
}

compare_betabinom <- function(coverage,ref_count,rho_val){
  sizes = coverage
  random = sapply(sizes, function(zz) rbetabinom(n=1, size=zz, prob = median_ratio ,rho =rho_val))
  
  suppressWarnings({ ks.test(random/sizes, ref_count/sizes)$statistic})
}

#===========================

adj_stat_all[1:3,]
row_id=which(gtex_w_ov$SAMPLE_ID=="GTEX-NPJ8-2626-SM-2D7W2")

ss<-gtex_w_ov %>% mutate(cut=cut_number(overlap,7)) %>% group_by(cut) %>% sample_n(5)
                       
#gtex_w_ov$SAMPLE_ID[which(gtex_w_ov$overlap>50)]
to_run<-ss$SAMPLE_ID
final<-data.frame(sample_id=rep(NA,length(to_run)), 
                  overlap=rep(NA,length(to_run)),
                  test_to_use=rep(NA,length(to_run)), 
                  disp_val=rep(NA,length(to_run)))

for(k in 1:length(to_run)){
  row_id=which(gtex_w_ov$SAMPLE_ID==to_run[k])
  
final$sample_id[k]=sample_id=gtex_w_ov$SAMPLE_ID[row_id]
study<-gtex_w_ov$study[row_id]
final$overlap[k]<-gtex_w_ov$overlap[row_id]
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
  }
}
    else{
      ase_df<-fread(saveTempFile)
    }
    #Fix the counts by adjusting the MA plot to mean around 0
    ratio_adj<-median(ase_df$ratio)
    #if(ratio_adj!=0){
      
      
      adj<-ase_df$ratio-ratio_adj
      ase_df$adj_alt<-round((ase_df$coverage/((2^adj)+1)),0)
      
      stopifnot(sum(is.na(ase_df$adj_alt))==0)
      
      ase_df$adj_ref<-ase_df$coverage-ase_df$adj_alt
    #}


#===========================
ase_df<-ase_df %>% filter(!adj_ref/coverage>0.85,!adj_ref/coverage<0.15)




#----------------------------------------------
median_ratio=median(ase_df$adj_ref/ase_df$coverage)
rho_seq<-seq(0,0.3,0.01)
betabinom_result<-data.frame(rho_seq)
betabinom_result$k.s.val<-NA
continue<-TRUE
for(i in 1:nrow(betabinom_result)){
  if(continue){
  print(i)
  out = replicate(10, compare_betabinom(ase_df$coverage,ase_df$adj_ref,rho_val=rho_seq[i]))
  random0 = sapply(ase_df$coverage, function(zz) rbetabinom(size = zz, n = 1, prob = median_ratio ,rho =rho_seq[i]))
  out0 = replicate(10, compare_betabinom(ase_df$coverage,random0,rho_val=rho_seq[i]))
  
  betabinom_result$k.s.val[i]<-median(out)-median(out0)
  
  if(sum(i>1, betabinom_result$k.s.val[i]>betabinom_result$k.s.val[i-1]) ==2){ continue<-FALSE}
 
  }
}
beta_binom_val<-min(betabinom_result$k.s.val,na.rm=T)
disper_val<-betabinom_result$rho_seq[which.min(betabinom_result$k.s.val)]

out = replicate(10, compare_binom(ase_df$coverage,ase_df$adj_ref))
random0 = sapply(ase_df$coverage, function(zz) rbinom(size = zz, n = 1, p = 0.5))
out0 = replicate(10, compare_binom(ase_df$coverage,random0))
binom_result<-median(out)-median(out0)


if(binom_result<beta_binom_val){  final$test_to_use[k]<-"binomial" 
}else{ final$test_to_use[k]<-"beta_binomial"
final$disp_val[k]<-disper_val}

}



test_ref = sapply(1:nrow(ase_df), function(zz) rbetabinom(n=1, size=ase_df$coverage[zz], prob = median_ratio,rho =0.02))
binom_ref = sapply(1:nrow(ase_df), function(zz) rbinom(n=1, size=ase_df$coverage[zz], prob = median_ratio ))
experimental_binom<- binom_ref/ase_df$coverage
experimental<- test_ref/ase_df$coverage
emperical<- ase_df$adj_ref/ase_df$coverage


df<-data.frame(experimental,emperical,experimental_binom)

pdf(file="~/plot/ASE/atest1.pdf", width = 10, height = 6)


ggplot(data=df ,aes(emperical))+
  geom_histogram(alpha=0.5)+
  geom_freqpoly(data=df,aes(experimental), color="darkgreen")+
  geom_freqpoly(data=df,aes(experimental_binom), color="darkred")+
  labs(title="Grey is emperical. Green is beta-binomial (disp=0.02), red is binomial.",
       subtitle= sample_id)
dev.off()

is.median<-median(ase_df$adj_ref/ase_df$coverage)
  ase_df$p_val = apply(ase_df[,c("adj_ref","adj_alt")], 1, function(x) {
    binom.test(round(x[1],1),round((x[1]+x[2]),1),p=is.median)$p.value})
  ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")
  
  ase_df$p_val_beta = apply(ase_df[,c("adj_ref","adj_alt")], 1, function(x) {
    2*(1-pbetabinom(x[1],(x[1]+x[2]),p=is.median,rho =0.02))})#rho_seq[i]
  id<-which(ase_df$adj_alt>ase_df$adj_ref)
  ase_df$p_val_beta[id] = sapply(id, function(xx) {
    pbetabinom(ase_df$adj_ref[xx],ase_df$coverage[xx],p=is.median,rho =0.02)*2 })
  
  ase_df$q_val_beta = p.adjust(ase_df$p_val_beta, method = "fdr")
  
  
  joined<-ase_df %>% 
    dplyr::select(c(chr, start,ref_count,alt_count,adj_alt, adj_ref,coverage,ref_count,alt_count,p_val,q_val)) %>% 
    inner_join(gtex_ase_1,by = c("chr", "start")) 
  
  length(unique(joined$GENE_ID[joined$q_val<0.05&joined$BINOM_P_ADJUSTED<0.05]))
  
  joined[1,]
sum(ase_df$q_val<0.05)
sum(gtex_ase_1$BINOM_P_ADJUSTED<0.05,na.rm=T)
sum(joined$q_val<0.05&joined$BINOM_P_ADJUSTED<0.05,na.rm=T)


pdf(file="~/plot/ASE/atest2.pdf", width = 10, height = 6)

ggplot(gtex_ase_1, aes(y=log2(ALT_COUNT),x=log2(REF_COUNT)))+
  geom_point(alpha=0.5)+
  geom_point(data=gtex_ase_1 %>% filter(BINOM_P_ADJUSTED<0.05), aes(y=log2(ALT_COUNT),x=log2(REF_COUNT)),color="red")

ggplot(ase_df, aes(y=log2(adj_alt),x=log2(adj_ref)))+
  geom_point(alpha=0.5)+
  geom_point(data=ase_df %>% filter(q_val_beta<0.05), aes(y=log2(adj_alt),x=log2(adj_ref)),color="red")

ggplot(ase_df, aes(y=log2(adj_alt),x=log2(adj_ref)))+
  geom_point(alpha=0.5)+
  geom_point(data=ase_df %>% filter(q_val<0.05), aes(y=log2(adj_alt),x=log2(adj_ref)),color="red")


ggplot(joined, aes(y=log2(adj_alt),x=log2(adj_ref)))+
  geom_point(alpha=0.5)+
  geom_point(data=joined %>% filter(q_val<0.05& BINOM_P_ADJUSTED<0.05), aes(y=log2(adj_alt),x=log2(adj_ref)),color="red")

dev.off()


xx<-ase_df %>% group_by(GENE_ID) %>% mutate(sig_cov=if_else(sum(q_val<0.05) >=1, min(coverage), 0)) %>% 
  summarize(n=n(),sig_percent=sum(q_val<0.05)/n, avg_cov=mean(coverage), min_sig_cov=unique(sig_cov)) %>% filter(sig_percent>0 ) %>% 
  arrange(desc(sig_percent))

ase_df %>% filter(GENE_ID=="ENSG00000002587") %>% tail(3)

xx %>% filter(n>1)%>% arrange(desc(sig))
