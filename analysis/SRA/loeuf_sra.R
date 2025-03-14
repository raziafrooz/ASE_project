library(data.table)
library(tidyverse)
library(scattermore)
library(ggridges)
source("~/ASE_project/src/remove_problematic_SNPs.R")
#source("~/ASE/src/uniNorm_function.R")

metadata<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/ASE_metadata.csv")
raw_metdata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/metadata/experiment_metadata.csv")
final_tissue<-fread("/dcs07/hansen/data/recount_ASE/data/toUse_tissue_primaryCell.csv")
uni_norm_df<-fread("/dcs07/hansen/data/recount_ASE/data/ase_uniNorm.csv")

ontology<-readRDS("/dcs07/hansen/data/recount_ASE/data/sra_ontology_term.rds")
prashanthi_annotations<-readRDS("~/ASE-data/data/prashanthi_annotations.rds")

sample_plot<-uni_norm_df[which(uni_norm_df$uni_norm_mean <= 0.085),]
final_tissue<- final_tissue %>% filter(experiment_acc %in% sample_plot$experiment_acc ) 
final_tissue<-final_tissue[-which(final_tissue$library_layout=="single" & final_tissue$bc_frag.mode_length!=0),]
#--------------------------------------------------

healthy<-final_tissue %>% filter(disease=="")

all_ase<-c()
for(i in 1454:1500){
  print(i)
exp_id<-healthy$experiment_acc[i]
ase_df<-fread(metadata$ASE_path[which(metadata$experiment_acc==exp_id)])

ase_df<-ase_df%>% 
  filter(pred_genotype==2, coverage>=8) %>%
  dplyr::select(chr,start,coverage,ref_count,alt_count,ref_ratio,pred_genotype,predicted_accuracy) %>%  
  mutate(ratio=log2(ref_count/alt_count),
         mean=(log2(ref_count)+log2(alt_count))/2)

#21047
ase_df<-remove_problematic_SNPs(ase_df=ase_df)
ase_df$genotyping_conf<-sapply(1:nrow(ase_df), function(zz) get_error_p(ase_df$ref_count[zz],ase_df$coverage[zz]))
ase_df<-ase_df %>% rowwise() %>% mutate(min_allele=min(ref_count,alt_count)) %>% filter(min_allele>2) %>% dplyr::select(!min_allele)

#ase_df_org<-ase_df
ase_df<-ase_df[which(ase_df$genotyping_conf>=0.000001),]


#Fix the counts by adjusting the MA plot to mean around 0 
ratio_adj<-median(ase_df$ratio)
adj<-ase_df$ratio-ratio_adj
ase_df$adj_alt<-round((ase_df$coverage/((2^adj)+1)),0)

stopifnot(sum(is.na(ase_df$adj_alt))==0)

ase_df$adj_ref<-ase_df$coverage-ase_df$adj_alt



is.median<-median(ase_df$adj_ref/ase_df$coverage)
ase_df$p_val = apply(ase_df[,c("adj_ref","adj_alt")], 1, function(x) {
  binom.test(round(x[1],1),round((x[1]+x[2]),1),p=is.median)$p.value})

ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")
ase_df$exp_id<-exp_id




ase_df<-ase_df %>% group_by(exp_id,symbol,chr,start) %>% summarize(n_sig_snps=sum(q_val<0.05))
all_ase<-rbind(all_ase,ase_df)

}

#fwrite(all_ase,"~/test/all_ase.csv")
#fwrite(all_ase,"~/test/all_ase2.csv")
#all_ase<-fread("~/test/all_ase.csv")
all_ase2<-fread("~/test/all_ase2.csv")
all_ase<-rbind(all_ase,all_ase2)

pdf(file="~/plot/ASE/sra.pdf", width = 10, height = 6)
ggplot(data=ase_df ,aes(y=log2(adj_alt), x=log2(adj_ref)))+
  geom_point(alpha=0.8)+
  geom_abline(slope=1, color="red")+
  geom_point(data=ase_df %>% filter(genotyping_conf<0.00001),aes(y=log2(adj_alt), x=log2(adj_ref)),alpha=0.8, color="red")+
  labs(title="233 samples in 4 tissues. blue is homo alt, purple is homo ref")

ggplot(data=ase_df ,aes(y=log2(adj_alt), x=log2(adj_ref)))+
  geom_point(alpha=0.8)+
  geom_abline(slope=1, color="red")+
  geom_point(data=ase_df %>% filter(q_val <0.05),aes(y=log2(adj_alt), x=log2(adj_ref)),alpha=0.8, color="red")+
  labs(title="233 samples in 4 tissues. blue is homo alt, purple is homo ref")

dev.off()
#=======================================================
gnomad<-fread("~/ASE-data/data/gnomad.v4.1.constraint_metrics.tsv")

gnomad<-gnomad %>%filter(mane_select==TRUE) %>%  dplyr::select(gene,lof.oe_ci.lower)
gnomad<-gnomad[!duplicated(gnomad$gene),]

total_samp<-length(unique(all_ase$exp_id))
all_ase_plot<-all_ase %>% group_by(symbol,exp_id) %>% 
  summarize(total_snp=n(),sig_snp=sum(n_sig_snps>0)) 

all_ase_plot2<-all_ase_plot %>% group_by(symbol) %>% 
  summarize(n_people=sum(sig_snp>0),perc=(sum(sig_snp>0)/total_samp)*100)%>% 
  arrange(desc(perc)) #%>% mutate(variability=cut_number(perc,3))

colnames(all_ase_plot2)[1]<-"gene"

all_ase_plot2<-left_join(all_ase_plot2,gnomad,by="gene")

all_ase_plot2<-all_ase_plot2 %>% mutate(variability=cut_number(perc,4),LOEUF=cut_number(lof.oe_ci.lower,10))
library(ggridges)
pdf(file="~/plot/ASE/sra_test2.pdf", width = 10, height = 6)
ggplot(all_ase_plot2, aes(x = lof.oe_ci.lower, y =variability,fill=variability)) + 
  geom_density_ridges2(alpha=0.9)+
  geom_vline(xintercept=0.35)

ggplot(all_ase_plot2, aes(x = perc, y =LOEUF,fill=LOEUF)) + 
  geom_density_ridges2(alpha=0.9)#+

dev.off()

all_ase_plot %>% filter(perc>0)
all_ase_plot$variability<-"not_ASE"
all_ase_plot$variability[all_ase_plot$perc>0.05]<-"ASE"

xx<-all_ase_plot %>% filter(variability=="ASE") %>% mutate(mm=cut_number(perc,7))


pdf(file="~/plot/ASE/sra_test2.pdf", width = 10, height = 6)
ggplot(xx, aes(x = lof.oe_ci.lower, y =mm,fill=mm)) + geom_density_ridges2()
dev.off()


all_ase_plot<-all_ase %>% group_by(symbol) %>% summarize(perc=(sum(n_sig_snps>0)/total_samp)*100)%>% 
  arrange(desc(perc)) %>%
  mutate(perc_cat=case_when(
    perc>3 ~ "high",
    perc<=3 ~ "low"
  )) %>% 
  filter(perc_cat=="high") %>% 
  mutate(variability=cut_number(perc,10))

colnames(all_ase_plot)[1]<-"gene"

all_ase_plot<-left_join(all_ase_plot,gnomad,by="gene")

all_ase_plot<-all_ase_plot %>% mutate(LOEUF=cut_number(lof.oe_ci.lower,10))

library(biomaRt)
snpmart = useEnsembl(biomart = "snp", dataset="hsapiens_snp")
listAttributes(snpmart)[1:10,1:2]
#clinical_significance
# rsid=c("rs123","rs150")
# 
# x<-ase_df[1:10,]
# xx<-getBM(attributes = c("chr_name","chrom_start","chrom_end","refsnp_id","allele",'consequence_type_tv'), 
#       filters = c('chr_name','start','end'), 
#       values = list(rep(1,2),x$start, x$start), 
#       mart = snpmart)
all_ase[duplicated(all_ase$chr & all_ase$start),]
all_ase[1:3,]
all_ase$loc<-paste0(all_ase$chr,"_", all_ase$start)
snps<-all_ase[!duplicated(all_ase$loc),] 

snps<-snps %>% arrange(chr)
chr_names<-str_replace(snps$chr,"chr","")

snpdb<-c()
for(ii in 1:40){
  print(ii)
  chr_name<-chr_names[ii]
  xx<-getBM(attributes = c("chr_name","chrom_start","chrom_end","refsnp_id","allele",'consequence_type_tv'), 
            filters = c('chr_name'), 
            values = list(chr_name), 
            mart = snpmart)
  snpdb<-rbind(snpdb,xx)
}
snpdb_filter<-snpdb %>% filter(chrom_start%in%snps$start)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

listAttributes(ensembl)[70:80,1:2]

xx<-getBM(attributes = c("start_position",'gene_biotype',"phenotype_description"), 
          filters = c('chromosome_name','start','end'), 
          values = list(8, 148350, 148420), 
          mart = ensembl)


#---------------------------------------------------------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
exon_gr<-exons(TxDb.Hsapiens.UCSC.hg38.knownGene)
t_gr<-transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

gnomad<-gnomad %>%filter(canonical==TRUE) %>%  dplyr::select(gene_id,gene,lof.oe_ci.upper,transcript_type,num_coding_exons)
gnomad<-gnomad[!is.na(gnomad$transcript_type),]
colnames(gnomad)[2]<-"symbol"
#length(unique(all_ase$exp_id)) #1453
all_ase_gr<-makeGRangesFromDataFrame(all_ase,seqnames="chr",start.field ="start",end.field = "start")
ov<-findOverlaps(all_ase_gr,exon_gr)
all_ase$IsExon<-NA
all_ase$IsExon[queryHits(ov)]<-TRUE
all_ase_exon<-all_ase %>% filter(IsExon==TRUE)
# all_ase_gr<-makeGRangesFromDataFrame(all_ase_exon,seqnames="chr",start.field ="start",end.field = "start")
# ov<-findOverlaps(all_ase_gr,t_gr)


n_snps_inExon<-all_ase_exon %>% group_by(symbol,exp_id) %>%  summarize(n_snps=n()) %>% ungroup()
n_snps_inExon<-n_snps_inExon %>% group_by(symbol) %>%  summarize(max_snps_inExon=max(n_snps),n_het_indv_inGene=length(unique(exp_id)))

n_snps_inExon %>% arrange(desc(max_snps_inExon)) # do we have more snps in genes that are larger?






summary_ase<-all_ase_exon %>% group_by(chr,start ,symbol) %>%
  summarize(n_ase=sum(n_sig_snps>0,na.rm=T),n_non_ase= sum(n_sig_snps==0,na.rm=T))
summary_ase<-left_join(summary_ase,n_snps_inExon,by=c("symbol")) %>% mutate(total_indv_inSNP=n_ase+n_non_ase)


summary_ase %>% arrange(desc(total_indv_inSNP))
summary_ase %>% filter(symbol=="RPS16") %>% arrange(desc(total_indv_inSNP))

# test<-summary_ase %>% group_by(symbol) %>% arrange(desc(n_ase)) %>% slice_head(n=1) %>% ungroup()
# test<-left_join(test,gnomad,by="symbol")

m<-summary_ase %>% filter(n_het_indv_inGene<10)
length(unique(m$symbol))
summary_ase %>% filter(symbol=="TMEM88B")
# 
xx2<-summary_ase %>% filter(n_het_indv_inGene>10) %>% group_by(symbol) %>% summarize(total_ase=sum(n_ase), total_non_ase=sum(n_non_ase))
xx2<-left_join(xx2,gnomad,by="symbol")
# 
# 
 xx2$perc_ase<-xx2$total_ase/(xx2$total_non_ase+xx2$total_ase)
 xx2$perc_non_ase<-xx2$total_non_ase/(xx2$total_non_ase+xx2$total_ase)
# xx2$total_snp<-(xx2$total_ase+xx2$total_non_ase)



never<-xx2 %>% filter(total_ase==0)
dim(never)

never<-xx2 %>% filter(perc_non_ase>0.98)
dim(never)


xx2<-xx2 %>%
  mutate(sig_ase_cut=cut_number(perc_ase,4),
         Loeuf_cut=cut_number(lof.oe_ci.upper,5))

pdf(file="~/plot/ASE/sra_test4.pdf", width = 8, height = 4)
ggplot(never, aes(lof.oe_ci.upper)) +
  geom_histogram()+labs(title="LOEUF for SNPs that never have sig ase")+
  geom_vline(xintercept=0.6,color="red")
ggplot(xx2, aes(y=sig_ase_cut,x=lof.oe_ci.upper,fill=sig_ase_cut)) +
  geom_density_ridges2()+
  geom_vline(xintercept=0.6)+
  labs(title="percent sig ase snps vs loeuf")+
  geom_vline(xintercept=0.6,color="red")

ggplot(xx2, aes(x=perc_ase,y=Loeuf_cut)) +
  geom_density_ridges2()+xlim(c(0,0.25))
dev.off()

#================================

Sys.setenv(ANNOTATION_HUB_CACHE=path.expand(rappdirs::user_cache_dir(appname="AnnotationHub")))
#get the gene locations that are protein coding only
ah <- AnnotationHub()
info <- query(ah, c("Homo.sapiens","Ensembl","GRCh38","gtf"))
grch38 <- info[['AH51014']]
#save(grch38, file = "/users/swang1/geuvadis_count/gene_anno_hg38.rda")
# subset for genes and chr1
chr_info <- c(1:22)
gene_info_hg38 <- grch38[!grch38$type %in% c('gene',"transcript")]
gene_info_hg38 <-gene_info_hg38[,1:10]
seqlevels(gene_info_hg38)<-paste0("chr",seqlevels(gene_info_hg38))

all_ase_gr<-makeGRangesFromDataFrame(summary_ase,seqnames="chr",start.field ="start",end.field = "start")
ov<-findOverlaps(all_ase_gr,gene_info_hg38)

summary_ase$type<-NA
summary_ase$type[queryHits(ov)]<-as.character(gene_info_hg38$type[subjectHits(ov)])

data.frame(summary_ase)[1:12,]
gene_info_hg38[1,]


xx2<-summary_ase %>% filter(n_het_indv_inGene>10) %>% group_by(symbol,type) %>% 
  summarize(total_ase=sum(n_ase), total_non_ase=sum(n_non_ase)) %>% ungroup()
xx2<-left_join(xx2,gnomad,by="symbol")
# 
# 
xx2$perc_ase<-xx2$total_ase/(xx2$total_non_ase+xx2$total_ase)
xx2$perc_non_ase<-xx2$total_non_ase/(xx2$total_non_ase+xx2$total_ase)
# xx2$total_snp<-(xx2$total_ase+xx2$total_non_ase)



never<-xx2 %>% filter(total_ase==0)
dim(never)


xx2<-xx2 %>%
  mutate(Loeuf_cut=cut_number(lof.oe_ci.upper,5))

pdf(file="~/plot/ASE/sra_test5.pdf", width = 8, height = 4)
ggplot(never, aes(lof.oe_ci.upper,fill=type)) +
  geom_histogram(alpha=0.8)+labs(title="LOEUF for SNPs that never have sig ase")+
  geom_vline(xintercept=0.6,color="red")


ggplot(xx2 %>% filter(!is.na(Loeuf_cut),!is.na(type) ), aes(x=total_ase,y=Loeuf_cut,fill=type)) +
  geom_density_ridges2(alpha=0.7)+xlim(c(0,100))

ggplot(xx2, aes(x=lof.oe_ci.upper,y=type,fill=type)) +
  geom_density_ridges2(alpha=0.7)+
  geom_vline(xintercept=0.6,color="red")
dev.off()



