library(data.table)
library(tidyverse)
source("~/ASE/src/uniNorm_function.R")

metadata<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/ASE_metadata.csv")
test_line<-readRDS("~/ASE/geuvadis_quantile_new.rds")
final_tissue<-fread("/dcs07/hansen/data/recount_ASE/data/toUse_tissue_primaryCell.csv")
uni_norm_df<-fread("/dcs07/hansen/data/recount_ASE/data/ase_uniNorm.csv")
ontology<-readRDS("/dcs07/hansen/data/recount_ASE/data/sra_ontology_term.rds")
prashanthi_annotations<-readRDS("~/ASE-data/data/prashanthi_annotations.rds")

sample_plot<-uni_norm_df[which(uni_norm_df$uni_norm_mean <= 0.085),]
final_tissue<- final_tissue %>% filter(experiment_acc %in% sample_plot$experiment_acc ) 
final_tissue<-final_tissue[-which(final_tissue$library_layout=="single" & final_tissue$bc_frag.mode_length!=0),]
#--------------------------
#no annotation:
#--------------------------
final_no_annot<-fread("/dcs07/hansen/data/recount_ASE/data/toUse_no_annot.csv")
uni_norm_df_noann<-fread("/dcs07/hansen/data/recount_ASE/data/ase_uniNorm_noAnno.csv")

sample_plot<-uni_norm_df_noann[which(uni_norm_df_noann$uni_norm_mean <= 0.085),]
final_no_annot<- final_no_annot %>% filter(experiment_acc %in% sample_plot$experiment_acc ) 
final_no_annot<-final_no_annot[-which(final_no_annot$library_layout=="single" & final_no_annot$bc_frag.mode_length!=0),]

final_no_annot$sample_type<-prashanthi_annotations$tissue.category[match(final_no_annot$experiment_acc,prashanthi_annotations$sample)]
final_no_annot$cancer<-prashanthi_annotations$cancer[match(final_no_annot$experiment_acc,prashanthi_annotations$sample)]
final_no_annot$cancer[is.na(final_no_annot$cancer)]<-"unknown"
final_no_annot<-final_no_annot %>% filter(cancer!="cancer",
                                          !sample_type %in% c("adipose tissue–derived MSCs",
                                                              "cancer cell line",
                                                              "cancer cell lines",
                                                              "Cells_Leukemia_cell_line_CML",
                                                              "day 10 hESCs cardiac lineage",
                                                              "day 6 hESCs cardiac lineage",
                                                              "day 2 hESCs cardiac lineage",
                                                              "differentiated iPSCs",
                                                              "eGFP- cells from Kidney organoid",
                                                              "eGFP+ cells from Kidney organoid",
                                                              "embryonic carcinoma",
                                                              "embryonic kidney",
                                                              "Hep2 cell line",
                                                              "hESC",
                                                              "hESC-derived cells",
                                                              "hESC differentiating into cardiomyocytes",
                                                              "hESC-derived cells",
                                                              "hESCs",
                                                              "hESCs presomitic mesoderm",
                                                              "hESCs somite",
                                                              "hiPSC",
                                                              "hiPSC-derived cells",
                                                              "iPSC",
                                                              "ipSCs",
                                                              "IPSCs",
                                                              "iPSCs",
                                                              "leukemia cell line",
                                                              "NK cells",
                                                              "reprogramming cells",
                                                              "reprogramming intermediates of hiF cells",
                                                              "reprogramming intermediates of hiF-T cells",
                                                              "sFRP2 knockdown mesenchymal stem cells",
                                                              "sFRP2OE mesenchymal stem cells",
                                                              "sFRP2 knockdown osteoblasts",
                                                              "sFRP2OE osteoblasts",
                                                              "transformed mesenchymal stem cells",
                                                              "universal human reference",
                                                              "universal human reference RNA",
                                                              "umbilical vein endothelial cells"))


#--------------------------
xx<-readRDS("~/ASE/tissue.rds")
ontology$parent<-xx
tissue_UBE<-ontology[str_starts(ontology$term_id, "UBERON"),]

tiss_id<-"kidney"
tiss_name<-unique(tissue_UBE$sample_accession[tissue_UBE$parent==tiss_id])


tissue<-final_tissue %>% filter(sample_accession %in% tiss_name)


#---------------------------------------
#find disease
#---------------------------------------


dim(final_tissue)

final_tissue$dd<-"No"
final_tissue$dd[final_tissue$disease!=""]<-"Yes"
table(final_tissue$dd,final_tissue$sample_type)

final_tissue[final_tissue$dd=="Yes",][5:10,]

xx<-ontology %>% filter(sample_accession %in% final_tissue$sample_accession,
                        str_starts(term_id, "DOID"))

xx %>% group_by(term_id) %>% summarize(n=n()) %>%  arrange(desc(n))


xx<-split(xx, xx$term_id )
names(xx)[15]
autism<-xx[["DOID:0060041"]] #autism
autism<-xx[["DOID:5419"]] #schizophrenia

autism<-final_tissue %>% filter(sample_accession %in% autism$sample_accession)
#---------------------------------------
pdf(file="~/plot/ASE/test2.pdf", width = 10, height = 6)

for(i in 1:nrow(paired)){
  print(i)
  exp_id<-paired$experiment_acc[i]
  sample_id<-paired$sample_accession[i]
  xx<-make_MA(exp_id, test_line)
  
  uni_norm<- xx$test_line_sample %>%
    filter(q=="high",CI=="high-ci",max >=1,max <=2.5) %>% 
    mutate(diff= ratio_q-CI_val) %>% 
    summarize(uni_norm_mean=mean(diff,na.rm=T))
  
  
  #uni_norm_df<-rbind(uni_norm_df,uni_norm)
  
  pp<-plot_MA(xx$ase_df, exp_id, sample_id,uni_norm$uni_norm_mean, xx$test_line_sample, xx$q_line)
  
  print(pp)
}
dev.off()


#---------------------------
sig_all<-c()
for(i in 1:nrow(autism)){
  print(i)
  exp_id<-autism$experiment_acc[i]
  ase_df<-fread(metadata$ASE_path[metadata$experiment_acc == exp_id])
  
  ase_df<-ase_df %>% 
    filter(coverage>=8, pred_genotype==2) %>% 
    mutate(alt_count=coverage-ref_count,
           ref_ratio=ref_count/coverage,
           log2aFC= log2((alt_count+1)/(ref_count+1)))
  
  median(ase_df$ref_ratio)
  
  
  
  #Make Granges:
  #-------------------------------
  ase_df_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")
  
  tx_38<-exons(TxDb,columns=c("TXNAME","GENEID","EXONID","CDSID"))
  
  ov<-findOverlaps(tx_38,ase_df_gr)
  
  ase_df$GENEID<-NA
  ase_df$GENEID[subjectHits(ov)]<-unlist(tx_38$GENEID[queryHits(ov)])
  
  
  
  ase_df<-ase_df %>% 
    dplyr::select(-c(AF,predicted_accuracy,pred_genotype,chr,start)) %>% 
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
    sig_all<-full_join(sig_all,sig, by="GENEID")
  }
  
}
#fwrite(sig_all,"~/test/autism.csv")
yy2<-ase_df %>% filter(GENEID=="ENSG00000184009.9")
library(matrixStats)
xx<-sig_all[,-1] %>% mutate(na=rowSums(is.na(.)),mm=rowMeans(., na.rm=T))

table(xx$na)
sig_all[which(xx$na<=350),][,c(1:3)]
xx$na[which(xx$na<25)]