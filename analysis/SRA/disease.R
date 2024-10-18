library(data.table)
library(tidyverse)
source("~/ASE/src/uniNorm_function.R")

metadata<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/ASE_metadata.csv")
final_tissue<-fread("/dcs07/hansen/data/recount_ASE/data/toUse_tissue_primaryCell.csv")
uni_norm_df<-fread("/dcs07/hansen/data/recount_ASE/data/ase_uniNorm.csv")
ontology<-readRDS("/dcs07/hansen/data/recount_ASE/data/sra_ontology_term.rds")
tx_38<-exons(TxDb,columns=c("TXNAME","GENEID","EXONID","CDSID"))

# Get LOEUF scores:
gnomad<-read.delim("hansen_lab/dwl_files/ASE_filter/gnomad.v4.1.constraint_metrics.tsv")

#GTEX ANEVA-dot Vg version 8
#load("/dcs07/hansen/data/recount_ASE/data/vg_ae.rda")
vg_ae<-read.delim("/dcs07/hansen/data/recount_ASE/data/Vg_GTEx_v8.txt")

sample_plot<-uni_norm_df[which(uni_norm_df$uni_norm_mean <= 0.085),]
final_tissue<- final_tissue %>% filter(experiment_acc %in% sample_plot$experiment_acc ) 


#---------------------------------------
#find disease
#---------------------------------------

xx<-ontology %>% filter(sample_accession %in% final_tissue$sample_accession,
                        str_starts(term_id, "DOID"))

xx %>% group_by(term_id) %>% summarize(n=n()) %>%  arrange(desc(n)) %>%  head(3)


xx<-split(xx, xx$term_id )

autism<-xx[["DOID:0060041"]] #autism
autism<-xx[["DOID:5419"]] #schizophrenia

autism<-final_tissue %>% filter(sample_accession %in% autism$sample_accession)
#---------------------------
i=1
sig_all<-c()
for(i in 42:nrow(autism)){
  print(i)
  exp_id<-autism$experiment_acc[i]
  ase_df<-fread(metadata$ASE_path[metadata$experiment_acc == exp_id])
  
  ase_df<-ase_df %>% 
    mutate(alt_count=coverage-ref_count,
           ref_ratio=ref_count/coverage)
  
  median(ase_df$ref_ratio)
  
  
  
  #Make Granges:
  #-------------------------------
  ase_df_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")
  
  ov<-findOverlaps(tx_38,ase_df_gr)
  
  ase_df$GENE_ID<-NA
  ase_df$GENE_ID[subjectHits(ov)]<-unlist(tx_38$GENEID[queryHits(ov)])
  
  ase_df$EXONID<-NA
  ase_df$EXONID[subjectHits(ov)]<-unlist(tx_38$EXONID[queryHits(ov)])
  # 
  
  ase_df$GENE_ID<-gsub("\\.\\d+", "", ase_df$GENE_ID)
  output_columns <- c("chr","start","GENE_ID","EXONID",  "REF_COUNT", "ALT_COUNT", "TOTAL_COUNT", "NULL_RATIO")
  #ase_df$TISSUE_ID<- "ADPSBQ"
  
  #xx<-ase_df %>% rowwise() %>% 
  #   mutate(hap_A=max(ref_count,alt_count), hap_B=min(ref_count,alt_count))
  
  test<-ase_df %>% 
    mutate(ALT_COUNT= coverage - ref_count, NULL_RATIO=median(ref_ratio),REF_COUNT = ref_count ,TOTAL_COUNT= coverage) %>% 
    dplyr::select(chr,start,GENE_ID,EXONID,REF_COUNT,ALT_COUNT, TOTAL_COUNT,NULL_RATIO ) %>% 
    filter(!is.na(GENE_ID))
  
  colnames(vg_ae)
  tiss <- "BRNCTXA" # The data comes from a skeletal muscle sample
  vg_ae_df<-vg_ae[match(test$GENE_ID,vg_ae$IDs),]
  covered_gene_Vgs <- vg_ae_df[,tiss]
  
  covered_gene_ASE_data <- test[match(vg_ae_df$IDs, test$GENE_ID),]
  
  
  # Take the square root of the Vg scores to the get the Standard Deviation (SDg)
  covered_gene_SDgs <- sqrt(covered_gene_Vgs) 
  
  
  # Run ANEVA-DOT
  ANEVADOT_recount <- ANEVADOT_test(covered_gene_ASE_data, output_columns = output_columns, 
                                    eh1 = "REF_COUNT", eh2 = "ALT_COUNT", coverage = 7, 
                                    r0 = covered_gene_ASE_data$NULL_RATIO,
                                    Eg_std = covered_gene_SDgs, plot = FALSE)
  
  ANEVADOT_recount$sample_id<-exp_id
  sig_all<-rbind(sig_all,ANEVADOT_recount)
}
write.csv(sig_all, "~/test/schizophernia.csv")

sig_all$log2aFC<-NULL
sig_all<-sig_all %>% mutate(log2aFC= (log2(ALT_COUNT+1)/log2(REF_COUNT+1)))

sig_all[which(sig_all$adj.pval<0.05),] %>% group_by(GENE_ID,EXONID) %>% summarize(n=n()) %>%  arrange(desc(n))





