#wasp<-gtex_tissue
#no_wasp<-gtex_tissue
library(AnnotationHub)
library("org.Hs.eg.db")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

ah <- AnnotationHub()
ah <- query(ah, c("v26","GENCODE","Homo sapiens","GRch38")) #v26 was used in GTExV8
TxDb<-ah[["AH75155"]]
tx_38<-exons(TxDb,columns=c("TXNAME","GENEID","EXONID","CDSID"))

tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("~/ASE-data/data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

tissues_names$file_name<-paste0("/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/", tissues_names$file_name)

gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)

#-----------------------------------------------------------------
#-----------------------------------------------------------------
xx<-read.csv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/metadata/all_gtex_metadata.csv")
gtex_sim<-fread("/dcs07/hansen/data/recount_ASE/data/gtex_simulation.csv.gz")
gtex_sim_gr<-makeGRangesFromDataFrame(gtex_sim,seqnames="chr",start.field ="start",end.field = "start")

recount_sig_snps<-c()

for (k in 3:nrow(xx)){
study<-xx$study[k]
print(study)
wasp_1<- fread(tissues_names$file_name[tissues_names$full_name==study][1]) %>% 
  filter(LOW_MAPABILITY<1,MAPPING_BIAS_SIM<1,GENOTYPE_WARNING<1)
colnames(wasp_1)[1:2]<- c("chr", "start")


xx_one<-xx[xx$study==study,]

#for(ss in 3:length(unique(xx$sample_id))){
for(ss in 1:length(unique(wasp_1$SAMPLE_ID))){
  print(ss)

sam<-unique(wasp_1$SAMPLE_ID)[ss]
sam_id<-xx_one$sample_id_rep[xx_one$sample_id==sam]


wasp_1_sam<-wasp_1 %>% filter(SAMPLE_ID==sam) 
if(nrow(wasp_1_sam)>1){
  


ase_df<-fread(gtex_metadata$genotypedSamples[gtex_metadata$sample_id_rep==sam][1]) %>%
  filter(pred_genotype==2, coverage>=8) %>% 
  mutate(ref_ratio=ref_count/coverage,
         alt_count=coverage-ref_count)



# perform multiple testing correction with FDR
#ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")


bigWig_path<-xx_one$total[which(xx_one$sample_id_rep==sam_id)]
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
  dplyr::select(chr, start,ref_count,alt_count, coverage,bigwig_count,ref_ratio )


ase_df$err_per <- (ase_df$bigwig_count - ase_df$coverage)/ase_df$bigwig_count


aa<-1-pbinom((ase_df$alt_count-1), size=ase_df$coverage, prob=mean(ase_df$err_per,na.rm=T))
rr<-1-pbinom((ase_df$ref_count-1), size=ase_df$coverage, prob=mean(ase_df$err_per,na.rm=T))


ase_df$geno_err<-rr+aa


ase_df<-ase_df[-which(ase_df$err_per>=0.05),] 
ase_df<-ase_df[-which(ase_df$geno_err>=0.001),]




ase_df_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")

ov<-findOverlaps(ase_df_gr,gtex_sim_gr)
ase_df$LOW_MAPABILITY[queryHits(ov)]<-gtex_sim$LOW_MAPABILITY[subjectHits(ov)]
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




is.median<-median(ase_df$ref_ratio)
ase_df$p_val = apply(ase_df[,c("ref_count","alt_count")], 1, function(x) {
  binom.test(round(x[1],1),round((x[1]+x[2]),1),p=is.median)$p.value})

ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")


wasp_1_sam<-wasp_1 %>% filter(SAMPLE_ID==sam) 
dd<-full_join(wasp_1_sam,ase_df,by=c("chr","start"))





dd_union<-inner_join(wasp_1_sam,ase_df,by=c("chr","start"))

n_union<-nrow(dd_union)
n_recount<-sum(!is.na(dd$ref_count))
n_gtex<-sum(!is.na(dd$REF_COUNT))
#before union:
sig_recount.05=sum(dd$q_val<0.05, na.rm=T)
sig_gtex= sum(dd$BINOM_P_ADJUSTED<0.05,na.rm=T)
gene_rec<-sum(!is.na((unique(dd$GENE_ID.y[dd$q_val<0.05]))))
gene_wasp<-sum(!is.na(unique(dd$GENE_ID.x[dd$BINOM_P_ADJUSTED<0.05])))
gene_both<-sum(!is.na(unique(dd$GENE_ID.x[dd$BINOM_P_ADJUSTED<0.05 & dd$q_val<0.05])))


sig_recount.05_uni=sum(dd_union$q_val<0.05, na.rm=T)
sig_gtex_uni= sum(dd_union$BINOM_P_ADJUSTED<0.05,na.rm=T)
sig_both_uni= sum(dd_union$BINOM_P_ADJUSTED<0.05 & dd_union$q_val<0.05,na.rm=T)
uni_gene_rec<-sum(!is.na((unique(dd_union$GENE_ID.y[dd_union$q_val<0.05]))))
uni_gene_wasp<-sum(!is.na(unique(dd_union$GENE_ID.x[dd_union$BINOM_P_ADJUSTED<0.05])))

df_x<-data.frame(SAMPLE_ID=sam,
                 total_n_union=n_union,
                 total_n_recount=n_recount,
                 total_n_gtex=n_gtex,
           sig_recount.05,sig_gtex,
           sig_recount.05_uni,sig_gtex_uni,sig_both_uni,
           gene_rec,gene_wasp,
           uni_gene_rec,uni_gene_wasp, gene_both,study)

recount_sig_snps<-rbind(recount_sig_snps,df_x)
}}
}
fwrite(recount_sig_snps, "~/test/recount_sig_snps.csv.gz")
# uni_gene_rec<-unique(ase_df$GENE_ID[ase_df$q_val<0.01])
# uni_gene_wasp<-unique(wasp_1_sam$GENE_ID[wasp_1_sam$BINOM_P_ADJUSTED<0.05])
# uni_gene_both<-unique(dd$GENE_ID[dd$BINOM_P_ADJUSTED<0.05 & dd$q_val<0.01])
# 
# 
# df_x<-data.frame(sample_id=sam,
#            sig_recount.01=sum(ase_df$q_val<0.01),
#            sig_recount.05=sum(ase_df$q_val<0.05),
#            sig_gtex= sum(wasp_1_sam$BINOM_P_ADJUSTED<0.05,na.rm=T),
#            sig_both= sum(dd$BINOM_P_ADJUSTED<0.05 & dd$q_val<0.05,na.rm=T),
#            only_recount= sum(dd$BINOM_P_ADJUSTED>=0.05 & dd$q_val<0.05,na.rm=T),
#            genes_recount=length(uni_gene_rec),
#            genes_wap=length(uni_gene_wasp),
#            genes_both=sum(uni_gene_rec %in% uni_gene_wasp)
#            )

#recount_sig_snps<-rbind(recount_sig_snps,df_x)

#}
#fwrite(recount_sig_snps, "~/test/recount_sig_snps.csv.gz")
#recount_sig_snps<-fread("~/test/recount_sig_snps.csv.gz")
qc_df<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/gtex_qc_metadata.csv.gz")
qc_df<-qc_df %>% mutate(overlap= (star.average_input_read_length)-bc_frag.mode_length)


#colnames(recount_sig_snps)[1]<-"SAMPLE_ID"

qc_df$SAMPLE_ID<-str_sub(qc_df$external_id, end= -3)
qc_df<-qc_df %>% select(SAMPLE_ID,overlap )
try1<-right_join(qc_df,recount_sig_snps)

try1 %>% filter(overlap>4) %>%  head(3)

plot_df<-try1 %>% sample_n(10) %>% 
  pivot_longer(!c(SAMPLE_ID,overlap), names_to = "pipeline", values_to = "n_sig")

#pdf(file="~/plot/ASE/significant_snps_wasp.pdf", width = 10, height = 6)
pdf(file="~/plot/ASE/test.pdf", width = 10, height = 6)


plot_df_1<-plot_df %>% filter(pipeline%in%c("sig_recount.05","sig_gtex"))

ggplot(data=plot_df_1, aes(x=SAMPLE_ID, y=n_sig, fill=pipeline)) +
  geom_bar(stat="identity", position=position_dodge())+
  labs(y="# of sig SNPs",
       title="# of SNPs with fdr<0.05 before union")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))

plot_df_1<-try1 %>% mutate(overlap_group=cut_number(overlap,3))%>% 
  pivot_longer(!c(SAMPLE_ID,overlap,overlap_group,study), names_to = "pipeline", values_to = "n_sig")%>% filter(pipeline%in%c("sig_recount.05","sig_gtex"))
ggplot(data=plot_df_1, aes(x=overlap_group, y=n_sig, fill=pipeline)) +
  geom_boxplot()+
  labs(y="# of sig SNPs",
       title="# of SNPs with fdr<0.05 before union")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))


plot_df_1<-plot_df %>% filter(pipeline%in%c("sig_recount.05_uni", "sig_gtex_uni", "sig_both_uni"))

ggplot(data=plot_df_1, aes(x=SAMPLE_ID, y=n_sig, fill=pipeline)) +
  geom_bar(stat="identity", position=position_dodge())+
  labs(y="# of sig SNPs",
       title="# of SNPs with fdr<0.05 union")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))

plot_df_1<-try1 %>% mutate(overlap_group=cut_number(overlap,3))%>% 
  pivot_longer(!c(SAMPLE_ID,overlap,overlap_group,study), names_to = "pipeline", values_to = "n_sig")%>% filter(pipeline%in%c("sig_recount.05_uni", "sig_gtex_uni", "sig_both_uni"))
ggplot(data=plot_df_1, aes(x=overlap_group, y=n_sig, fill=pipeline)) +
  geom_boxplot()+
  labs(y="# of sig SNPs",
       title="# of SNPs with fdr<0.05 before union")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))


plot_df_1<-plot_df %>% filter(pipeline%in%c("gene_rec","gene_wasp"))


ggplot(data=plot_df_1, aes(x=SAMPLE_ID, y=n_sig, fill=pipeline)) +
  geom_bar(stat="identity", position=position_dodge())+
  labs(y="# of sig genes",
       title="# of sig genes (recount fdr 0.05) before union")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))

plot_df_1<-try1 %>% mutate(overlap_group=cut_number(overlap,3))%>% 
  pivot_longer(!c(SAMPLE_ID,overlap,overlap_group,study), names_to = "pipeline", values_to = "n_sig")%>% filter(pipeline%in%c("gene_rec","gene_wasp"))
ggplot(data=plot_df_1, aes(x=overlap_group, y=n_sig, fill=pipeline)) +
  geom_boxplot()+
  labs(y="# of sig SNPs",
       title="# of SNPs with fdr<0.05 before union")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))



plot_df_1<-plot_df %>% filter(pipeline%in%c("uni_gene_rec","uni_gene_wasp","gene_both"))


ggplot(data=plot_df_1, aes(x=SAMPLE_ID, y=n_sig, fill=pipeline)) +
  geom_bar(stat="identity", position=position_dodge())+
  labs(y="# of sig genes",
       title="# of sig genes (recount fdr 0.05) after union")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))

plot_df_1<-try1 %>% mutate(overlap_group=cut_number(overlap,3)) %>% 
  pivot_longer(!c(SAMPLE_ID,overlap,overlap_group,study), names_to = "pipeline", values_to = "n_sig")%>% filter(pipeline %in%c("uni_gene_rec","uni_gene_wasp","gene_both"))
ggplot(data=plot_df_1, aes(x=overlap_group, y=n_sig, fill=pipeline)) +
  geom_boxplot()+
  labs(y="# of sig SNPs",
       title="# of SNPs with fdr<0.05 before union")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))


dev.off()

try1<-try1%>%filter(sig_gtex>0)
try1$overlap_group<-cut_number(try1$overlap,5)

pdf(file="~/plot/ASE/significant_snps_wasp2.pdf", width = 10, height = 6)

ggplot(data=try1, aes(x=overlap_group, y=sig_both)) +
  geom_boxplot()+
  labs(y="# same sig snps",
       title= "Number of SNPs that are significant both in wasp and recount")

ggplot(data=try1, aes(x=overlap_group, y=genes_both)) +
  geom_boxplot()+
  labs(y="# same sig genes",
       title= "Number of genes that are significant both in wasp and recount")

dev.off()







#-----------------------
pdf(file="~/plot/ASE/overlap_sig.pdf", width = 10, height = 6)


plot_df_1<-try1 %>% mutate(overlap_group=cut_number(overlap,3))%>% 
  pivot_longer(!c(SAMPLE_ID,overlap,overlap_group,study), names_to = "pipeline", values_to = "n_sig")%>% filter(pipeline%in%c("sig_recount.05","sig_gtex"))
ggplot(data=plot_df_1, aes(x=overlap, y=n_sig, color=pipeline)) +
  geom_point(alpha=0.4)+
  labs(y="# of sig SNPs",
       title="# of SNPs with fdr<0.05 before union")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))

plot_df_1<-try1 %>% mutate(overlap_group=cut_number(overlap,3))%>% 
  pivot_longer(!c(SAMPLE_ID,overlap,overlap_group,study), names_to = "pipeline", values_to = "n_sig")%>% filter(pipeline%in%c("gene_rec","gene_wasp"))
ggplot(data=plot_df_1, aes(x=overlap, y=n_sig, color=pipeline)) +
  geom_point()+
  labs(y="# of sig SNPs",
       title="# of SNPs with fdr<0.05 before union")+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))

dev.off()





