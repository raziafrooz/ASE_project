library(data.table)
library(tidyverse)
library(ANEVADOT)
library(AnnotationHub)
setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(MetBrewer)

source("~/ASE/src/uniNorm_function.R")


load("/dcs07/hansen/data/recount_ASE/data/vg_ae.rda")

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

# 
# metadata<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/ASE_metadata.csv")
# test_line<-readRDS("~/ASE/geuvadis_quantile_new.rds")
# final_tissue<-fread("/dcs07/hansen/data/recount_ASE/data/toUse_tissue_primaryCell.csv")
# uni_norm_df<-fread("/dcs07/hansen/data/recount_ASE/data/ase_uniNorm.csv")
# ontology<-readRDS("/dcs07/hansen/data/recount_ASE/data/sra_ontology_term.rds")


#------------------------------------------------------
#Get the true ASE hits from gtex 
#------------------------------------------------------
#met<-read.csv("data/GTEx_metadata.csv")

tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

tissues_names$file_name<-paste0("/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/", tissues_names$file_name)


#-----------------------------------------------------
#Compare true GTEx to Reocunt3
#-----------------------------------------------------
#Get gtex genotype metadata:
gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")



# Define the output columns
output_columns <- c("GENE_ID", "TISSUE_ID",  "REF_COUNT", "ALT_COUNT", "TOTAL_COUNT", "NULL_RATIO")

ase_df<-fread(gtex_metadata$genotypedSamples[1]) %>%
  filter(pred_genotype==2, coverage>=8) %>% 
  mutate(ref_ratio=ref_count/coverage)


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
#---------------------------------------------------------------
#get gene id


ase_df_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")

tx_38<-exons(TxDb,columns=c("TXNAME","GENEID","EXONID","CDSID"))

ov<-findOverlaps(tx_38,ase_df_gr)

ase_df$GENE_ID<-NA
ase_df$GENE_ID[subjectHits(ov)]<-unlist(tx_38$GENEID[queryHits(ov)])


ase_df$GENE_ID<-gsub("\\.\\d+", "", ase_df$GENE_ID)
output_columns <- c("chr","start","GENE_ID", "TISSUE_ID",  "REF_COUNT", "ALT_COUNT", "TOTAL_COUNT", "NULL_RATIO")
ase_df$TISSUE_ID<- "ADPSBQ"

test<-ase_df %>% 
  mutate(ALT_COUNT= coverage - ref_count, NULL_RATIO=median(ref_ratio),REF_COUNT = ref_count ,TOTAL_COUNT= coverage) %>% 
  dplyr::select(chr,start,GENE_ID,TISSUE_ID,REF_COUNT,ALT_COUNT, TOTAL_COUNT,NULL_RATIO ) %>% 
  filter(!is.na(GENE_ID))


tiss <- "MEAN" # The data comes from a skeletal muscle sample
covered_genes <- intersect(vg_ae$IDs, test$GENE_ID)
covered_gene_Vgs <- vg_ae[match(covered_genes, vg_ae$IDs), tiss] 
covered_gene_ASE_data <- test %>% filter(GENE_ID %in% covered_genes)

# Take the square root of the Vg scores to the get the Standard Deviation (SDg)
covered_gene_SDgs <- sqrt(covered_gene_Vgs) 


# Run ANEVA-DOT
ANEVADOT_recount <- ANEVADOT_test(covered_gene_ASE_data, output_columns = output_columns, 
                                 eh1 = "REF_COUNT", eh2 = "ALT_COUNT", coverage = 7, 
                                 r0 = covered_gene_ASE_data$NULL_RATIO,
                                 Eg_std = covered_gene_SDgs, plot = FALSE)


# pdf(file="~/plot/ASE/test4.pdf", width = 10, height = 6)
# ggplot(ANEVADOT_scores, aes(x=log2(REF_COUNT), y=log2(ALT_COUNT)))+
#   geom_point(alpha=0.5)+
#   geom_point(data=ANEVADOT_scores[which(ANEVADOT_scores$adj.pval<0.05)],aes(x=log2(REF_COUNT), y=log2(ALT_COUNT)),color="red")+
#   labs(title="AR ran the first GTEx sample from ADPSBQ tissue for test")
# 
# dev.off()


# xx<-ANEVADOT_scores[which(ANEVADOT_scores$adj.pval<0.05)]
# dim(xx)
# ANEVADOT_scores[
#   (match(unique(xx$GENE_ID),unique(xx2$GENE_ID)))


#---------------------------
#true GTEx

# ase_df<-fread(gtex_metadata$genotypedSamples[1]) %>%
#   filter(pred_genotype==2, coverage>=8) %>% 
#   mutate(ref_ratio=ref_count/coverage)

gtex_tissue<- fread(tissues_names$file_name[tissues_names$abb==tiss])
colnames(gtex_tissue)[1:2]<- c("chr", "start")

gtex_df<-gtex_tissue %>% filter(SAMPLE_ID=="GTEX-1117F-0226-SM-5GZZ7")

ase_df_gr<-makeGRangesFromDataFrame(gtex_df,seqnames="chr",start.field ="start",end.field = "start")

tx_38<-exons(TxDb,columns=c("TXNAME","GENEID","EXONID","CDSID"))

ov<-findOverlaps(tx_38,ase_df_gr)

gtex_df$GENE_ID<-NA
gtex_df$GENE_ID[subjectHits(ov)]<-unlist(tx_38$GENEID[queryHits(ov)])


gtex_df$GENE_ID<-gsub("\\.\\d+", "", gtex_df$GENE_ID)
output_columns <- c("chr","start","GENE_ID", "TISSUE_ID",  "REF_COUNT", "ALT_COUNT", "TOTAL_COUNT", "NULL_RATIO")
gtex_df$TISSUE_ID<- "ADPSBQ"

test<-gtex_df %>% dplyr::select(chr,start, GENE_ID,TISSUE_ID,REF_COUNT,ALT_COUNT, TOTAL_COUNT,NULL_RATIO ) %>% 
  filter(!is.na(GENE_ID))


tiss <- "MEAN" # The data comes from a skeletal muscle sample
covered_genes <- intersect(vg_ae$IDs, test$GENE_ID)
covered_gene_Vgs <- vg_ae[match(covered_genes, vg_ae$IDs), tiss] 
covered_gene_ASE_data <- test %>% filter(GENE_ID %in% covered_genes)

# Take the square root of the Vg scores to the get the Standard Deviation (SDg)
covered_gene_SDgs <- sqrt(covered_gene_Vgs) 


# Run ANEVA-DOT
ANEVADOT_gtex <- ANEVADOT_test(covered_gene_ASE_data, output_columns = output_columns, 
                                 eh1 = "REF_COUNT", eh2 = "ALT_COUNT", coverage = 7, 
                                 r0 = covered_gene_ASE_data$NULL_RATIO,
                                 Eg_std = covered_gene_SDgs, plot = FALSE)


# 
# pdf(file="~/plot/ASE/test3.pdf", width = 10, height = 6)
# ggplot(ANEVADOT_scores, aes(x=log2(REF_COUNT), y=log2(ALT_COUNT)))+
#   geom_point(alpha=0.5)+
#   geom_point(data=ANEVADOT_scores[which(ANEVADOT_scores$adj.pval<0.05)],aes(x=log2(REF_COUNT), y=log2(ALT_COUNT)),color="red")+
#   labs(title="AR ran the first GTEx sample from ADPSBQ tissue for test (True GTEx)")
# 
# dev.off()
# 



#-------------------------------
output_columns <- c("chr","start","GENE_ID", "TISSUE_ID", 
                    "REF_COUNT_recount", "ALT_COUNT_recount", "TOTAL_COUNT_recount", "NULL_RATIO_recount",
                    "p.val_recount", "adj.pval_recount")
colnames(ANEVADOT_recount)<-output_columns

dim(ANEVADOT_recount)




both<-full_join(ANEVADOT_recount,ANEVADOT_gtex)
both<-both %>% mutate(aFC_reocunt= log2(ALT_COUNT_recount+1)/log2(REF_COUNT_recount+1),
                aFC_gtex= log2(ALT_COUNT+1)/log2(REF_COUNT+1))

xx_gtex<-both[which(both$adj.pval<0.05)]
xx_recount<-both[which(both$adj.pval_recount<0.05)]

fwrite(xx_recount,"~/ANEVA_recount_sig.csv")
fwrite(xx_gtex,"~/ANEVA_gtex_sig.csv")

union<-both[which(both$adj.pval<0.05 & both$adj.pval_recount<0.05)]
fwrite(union,"~/ANEVA_union_sig.csv")

