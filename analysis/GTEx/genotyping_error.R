library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library("org.Hs.eg.db")
library(VGAM)
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

gtex_predict<-read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/metadata/all_gtex_metadata.csv")
gtex_sim<-fread("/dcs07/hansen/data/recount_ASE/data/gtex_simulation.csv.gz")
gtex_sim_gr<-makeGRangesFromDataFrame(gtex_sim,seqnames="chr",start.field ="start",end.field = "start")

gtex_test<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/GTEx_testing.csv")
gtex_uni_norm<-fread("/dcs07/hansen/data/recount_ASE/data/gtex_uniNorm.csv")

bad_samples<-gtex_uni_norm$sample_id[which(gtex_uni_norm$uni_norm_mean>0.08)]

#============================================================
#Read in wasp indv:
path_id<-"/dcl01/hansen/data/arazi/ASE/dbGap/phe000039.v1.GTEx_v8_ASE_WASP.expression-matrixfmt-ase.c1/GTEx_Analysis_v8_ASE_WASP_counts_by_subject/"
tissues_names<-as.data.frame(list.files(path =path_id))
colnames(tissues_names)<-"file_name"
tissues_names$indv<-gsub("\\..*","",tissues_names$file_name)
tissues_names$file_name<-paste0(path_id, tissues_names$file_name)
#============================================================
tissues<-unique(gtex_test$tissue)


for(t_id in 4:7){
  tissue_name<-tissues[t_id]
  print(tissue_name)
  
  if(!file.exists(paste0("~/test/geno_error_",tissue_name,".csv.gz"))){
   
gtex_prediction<-fread(gtex_test$genotypedSamples[gtex_test$tissue==tissue_name])
samples<-unique(gtex_prediction$sample_id_rep)
rm_id<-which(str_sub(samples, end= -3) %in% bad_samples)
if(length(rm_id)>0){samples<-samples[-rm_id]}

print(paste0("number of sample:",length(samples)))
all_ase<-c()
for(row_id in 1:length(samples)){
  print(row_id)
  sample_id<-samples[row_id]
  
    ase_df<-gtex_prediction%>% 
      filter(sample_id_rep== sample_id,pred_genotype==2, coverage>=8) %>% 
      mutate(ref_ratio=ref_count/coverage,
             alt_count=coverage-ref_count,
             ratio=log2(ref_count/alt_count),
             mean=(log2(ref_count)+log2(alt_count))/2)
    
    
        bigWig_path<-gtex_predict$total[which(gtex_predict$sample_id_rep==sample_id)][1]
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
        
     
      #Fix the counts by adjusting the MA plot to mean around 0
      ratio_adj<-median(ase_df$ratio)
      #if(ratio_adj!=0){
        
        
        adj<-ase_df$ratio-ratio_adj
        ase_df$adj_alt<-round((ase_df$coverage/((2^adj)+1)),0)
        
        stopifnot(sum(is.na(ase_df$adj_alt))==0)
        
        ase_df$adj_ref<-ase_df$coverage-ase_df$adj_alt
        
        ase_df<-ase_df %>%
          dplyr::select(chr, start,adj_alt,adj_ref, coverage,true_genotype,pred_genotype,sample_id_rep)
        
        all_ase<-rbind(all_ase,ase_df)   
     
    }

fwrite(all_ase,paste0("~/test/geno_error_",tissue_name,".csv.gz"))
}
}
all_ase<-fread("~/test/geno_error.csv.gz")

for(t_id in 4:7){
  tissue_name<-tissues[t_id]
  print(tissue_name)
  
  ase<-fread(paste0("~/test/geno_error_",tissue_name,".csv.gz"))
  
  all_ase<-rbind(all_ase,ase)
}

pdf(file="~/plot/ASE/geno_error.pdf", width = 10, height = 6)
ggplot(data=all_ase ,aes(y=log2(adj_alt), x=log2(adj_ref)))+
  geom_scattermore(alpha=0.8,pointsize=2 )+
  geom_abline(slope=1, color="red")+
  geom_scattermore(data=all_ase %>% filter(true_genotype==1),aes(y=log2(adj_alt), x=log2(adj_ref)),aalpha=0.8,pointsize=2, color="purple")+
  geom_scattermore(data=all_ase %>% filter(true_genotype==3),aes(y=log2(adj_alt), x=log2(adj_ref)),alpha=0.8,pointsize=2, color="blue")+
  labs(title="233 samples in 4 tissues. blue is homo alt, purple is homo ref")
  
ggplot(data=all_ase %>% filter(true_genotype==2),aes(adj_ref/coverage))+
  geom_histogram(alpha=0.4)+
  geom_histogram(data=all_ase %>% filter(true_genotype==1),aes(adj_ref/coverage),fill="purple",alpha=0.4)+
  geom_histogram(data=all_ase %>% filter(true_genotype==3),aes(adj_ref/coverage),fill="blue",alpha=0.4)+
  labs(title="233 samples in 4 tissues. gray is true het/purple is true homo ref/ blue is true homo alt")

dev.off()



compare_betabinom <- function(coverage,ref_count,rho_val){
  sizes = coverage
  random = sapply(sizes, function(zz) rbetabinom(n=1, size=zz, prob = median_ratio ,rho =rho_val))
  
  suppressWarnings({ ks.test(random/sizes, ref_count/sizes)$statistic})
}

xx<-all_ase %>% filter(true_genotype==3)
yy<-all_ase %>% filter(true_genotype==1)



median_ratio=median(yy$adj_ref/yy$coverage)
rho_seq<-seq(0,0.3,0.01)
betabinom_result<-data.frame(rho_seq)
betabinom_result$k.s.val<-NA
continue<-TRUE
for(i in 1:nrow(betabinom_result)){
  #if(continue){
    print(i)
    out = replicate(10, compare_betabinom(yy$coverage,yy$adj_ref,rho_val=rho_seq[i]))
    random0 = sapply(yy$coverage, function(zz) rbetabinom(size = zz, n = 1, prob = median_ratio ,rho =rho_seq[i]))
    out0 = replicate(10, compare_betabinom(yy$coverage,random0,rho_val=rho_seq[i]))
    
    betabinom_result$k.s.val[i]<-median(out)-median(out0)
    print(betabinom_result$k.s.val[i])
    #if(sum(i>1, betabinom_result$k.s.val[i]>betabinom_result$k.s.val[i-1]) ==2){ continue<-FALSE}
    
  #}
}

beta_binom_val<-min(betabinom_result$k.s.val,na.rm=T)
disper_val<-betabinom_result$rho_seq[which.min(betabinom_result$k.s.val)]


test_ref = sapply(1:nrow(yy), function(zz) rbetabinom(n=1, size=yy$coverage[zz], prob = median_ratio,rho =disper_val))
binom_ref = sapply(1:nrow(yy), function(zz) rbinom(n=1, size=yy$coverage[zz],median_ratio))
experimental_binom<- binom_ref/yy$coverage
experimental<- test_ref/yy$coverage
emperical<- yy$adj_ref/yy$coverage


df<-data.frame(experimental,emperical,experimental_binom)

pdf(file="~/plot/ASE/geno_error_test2.pdf", width = 10, height = 6)


ggplot(data=df ,aes(emperical))+
  geom_histogram(alpha=0.5)+
  geom_freqpoly(data=df,aes(experimental), color="darkgreen")+
  geom_freqpoly(data=df,aes(experimental_binom), color="darkred")+
  labs(title="Grey is emperical. Green is beta-binomial (disp=0.07)")
dev.off()


try<-all_ase #%>% group_by(true_genotype) %>% sample_n(5)
aa<-sapply(1:nrow(try), function(zz)  { pbetabinom((try$adj_ref[zz]), size=try$coverage[zz], prob=0.8276,rho=0) })
try<- try %>% rowwise() %>% mutate(min_allele=min(adj_alt,adj_ref))

rr<-try[-which(aa>0.3),]
table(try$true_genotype)
table(rr$true_genotype)

library(scattermore)
pdf(file="~/plot/ASE/geno_error_result.pdf", width = 10, height = 6)

ggplot(rr,aes(y=min_allele,x=coverage))+
  geom_scattermore(pointsize=2)+
  geom_scattermore(data=try %>% filter(true_genotype==1),aes(y=min_allele,x=coverage), pointsize=2,color="red")+
  geom_scattermore(data=rr %>% filter(true_genotype==1),aes(y=min_allele,x=coverage), pointsize=2,color="purple")+
  xlim(c(0,1000))+
  ylim(c(0,500))+
  labs(title="233 samples in 4 tissues.purple is true homo ref not removed/ red homo ref removed (geno-error>0.3)")


ggplot(rr,aes(y=log2(adj_alt),x=log2(adj_ref)))+
  geom_scattermore(pointsize=2)+
  geom_abline(slope=1, color="red")+
  geom_scattermore(data=rr %>% filter(true_genotype==1),aes(y=log2(adj_alt), x=log2(adj_ref)),pointsize=2, color="purple")+
  labs(title="233 samples in 4 tissues.purple is true homo ref not removed/ red homo ref removed (geno-error>0.)")

dev.off()
  


#================================
#ROC curve
vv=0.1
try$error_perc<-aa
rr<-try

roc_df<-data.frame(cut_off=c(0.01,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9))
for(i in 1:nrow(roc_df)){
  vv<-roc_df$cut_off[i]
  print(vv)
rr$rmv<-"no"
rr$rmv[rr$error_perc>=vv]<-"yes"


roc_df$true_pos[i]<-sum(rr$true_genotype==1 & rr$rmv=="yes") /sum(rr$true_genotype==1)
roc_df$false_pos[i]<-sum(rr$true_genotype==2 & rr$rmv=="yes") /sum(rr$true_genotype==2)

}

pdf(file="~/plot/ASE/geno_error_roc.pdf", width = 10, height = 6)
ggplot(roc_df, aes(x=false_pos,y=true_pos,group=1,label=cut_off))+
  geom_line()+
  geom_point()+geom_text(hjust=0, vjust=0)+
  labs(title="ROC curve for genotyping error cut-off")
dev.off()
