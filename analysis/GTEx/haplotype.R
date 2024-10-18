library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(MetBrewer)
library(AnnotationHub)
ah <- AnnotationHub()
ah <- query(ah, c("v26","GENCODE","Homo sapiens","GRch38")) #v26 was used in GTExV8
TxDb<-ah[["AH75155"]]

tx_38<-exons(TxDb,columns=c("TXNAME","GENEID","EXONID","CDSID"))


#Downloaded blacklist of snps from gtex (https://github.com/secastel/phaser/blob/master/phaser/README.md)
HLA_snp<-read.table("/users/arazi/hansen_lab/dwl_files/ASE_filter/hg38_haplo_count_blacklist.chr.bed", sep="\t")
HLA_snp_gr<-makeGRangesFromDataFrame(HLA_snp,seqnames="V1",start.field ="V2",end.field = "V3")

#I found a bed file including the problematic sites from ENCODE
#https://github.com/Boyle-Lab/Blacklist/tree/master
black_list<-import("/users/arazi/hansen_lab/dwl_files/ASE_filter/hg38-blacklist.v2.bed.gz")

#https://bismap.hoffmanlab.org
mappability<-import("/users/arazi/hansen_lab/dwl_files/ASE_filter/k50.Unique.Mappability.bb")


phaser<-fread("~/test/phASER_GTEx_v8_matrix_WASP.gw_phased.txt")
colnames(phaser)[1:2]<-c("chr","GENE_ID")

gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)



colnames(phaser)<-gsub("[.]","-",colnames(phaser))


for(i in 1:3){
id<-which(gtex_metadata$sample_id_rep %in% colnames(phaser))[i]

subject_id<-colnames(phaser)[id]
sam<-colnames(phaser)[id]



ase_df<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==sam][1]) %>%
  filter(pred_genotype==2, coverage>=8) %>% 
  mutate(ref_ratio=ref_count/coverage,
         alt_count= coverage - ref_count)

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


ov<-findOverlaps(tx_38,ase_filt_gr)

ase_df$GENE_ID<-NA
ase_df$GENE_ID[subjectHits(ov)]<-unlist(tx_38$GENEID[queryHits(ov)])


one_sample<-phaser[,c(1:4, id)]
colnames(one_sample)[5]<-"phaser"

jj<-left_join(ase_df,one_sample, by=c("chr","GENE_ID"))





numbers<-jj %>% group_by(GENE_ID) %>%  summarize(n=n()) %>% arrange(desc(n))
numbers<-numbers[numbers$n>1 & !is.na(numbers$GENE_ID),]

all_genes<-c()
for(k in 1:nrow(numbers)){

  print(k)
  g_id<-numbers$GENE_ID[k]
  
xx<-jj %>%
  filter(GENE_ID==g_id) %>%
  rowwise() %>% 
  mutate(hap_A=max(ref_count,alt_count), hap_B=min(ref_count,alt_count))

xx$hap_A= sum(xx$hap_A)
xx$hap_B=sum(xx$hap_B)

xx<-xx %>% dplyr::select(GENE_ID,hap_A, hap_B, phaser) %>% head(1)

all_genes<-rbind(all_genes,xx)

}



xx<-strsplit(all_genes$phaser,"\\|")
xx<-do.call(rbind.data.frame, xx)
colnames(xx)<-c("A", "B")
xx$A<-as.numeric(xx$A)
xx$B<-as.numeric(xx$B)
xx<-xx %>% rowwise() %>% mutate(hap_A=max(A,B), hap_B=min(A,B))


all_genes$phaser_hap_A<-as.integer(xx$hap_A)
all_genes$phaser_hap_B<-as.integer(xx$hap_B)

all_genes<-all_genes %>% rowwise() %>%
  mutate(ratioA=hap_A/(hap_A+hap_B),ratioA_phaser=phaser_hap_A/(phaser_hap_A+phaser_hap_B),
  ratio_hapA=log2(ratioA_phaser) - log2(ratioA),
                     mean_hapA=(log2(ratioA_phaser) + log2(ratioA))/2)#,
                    # ratio_hapB=log2(phaser_hap_B) - log2(hap_B),
                    # mean_hapB=(log2(phaser_hap_B) + log2(hap_B))/2)

pdf(file=paste0("~/plot/ASE/haplotype",subject_id,"_ratio2.pdf"), width = 10, height = 6)


pp=ggplot(all_genes, aes(y=ratio_hapA, x=mean_hapA))+
  geom_point(alpha=0.5)+
  geom_abline(slope=0, color="red")+
  labs(title=paste0(subject_id,", subseted for the SNPs that passes WASP on GTEx"))

print(pp)

pp=ggplot(all_genes, aes(y=ratio_hapB, x=mean_hapB))+
  geom_point(alpha=0.5)+
  geom_abline(slope=0, color="red")+
  labs(title=subject_id)

print(pp)

dev.off()

}

subject_id
path<-"/dcl01/hansen/data/arazi/ASE/dbGap/phe000039.v1.GTEx_v8_ASE_WASP.expression-matrixfmt-ase.c1/GTEx_Analysis_v8_ASE_WASP_counts_by_subject/"
wasp<-fread(paste0(path,"GTEX-12WSN.v8.wasp_corrected.ase_table.tsv.gz"))
wasp[1,1:4]
colnames(wasp)[1:2]<-c("chr", "start")
jj<-inner_join(ase_df, wasp, by=c("chr", "start"))

colnames(jj)[12]<-"GENE_ID"
jj<-jj[,1:ncol(ase_df)]


all_genes %>% filter(ratio_hapA>0.8, mean_hapA> -0.5)




