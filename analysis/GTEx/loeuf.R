library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(MetBrewer)    
library("AnnotationDbi")
library("org.Hs.eg.db")
library(ggridges)

tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("~/ASE-data/data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

tissues_names$file_name<-paste0("/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/", tissues_names$file_name)


#-----------------------------------------------------
#Compare true GTEx to Reocunt3
#-----------------------------------------------------
#Get gtex genotype metadata:
gtex_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/AggregateFiles/all_GTEx.csv")
gtex_metadata$sample_id_rep<-str_sub(gtex_metadata$sample_id, end= -3)

gnomad<-fread("~/ASE-data/data/gnomad.v4.1.constraint_metrics.tsv")

tissues<-tissues_names$full_name


  study<-tissues[1]
  print(study)
  

gtex_tissue<- fread(tissues_names$file_name[tissues_names$full_name==study][1])
gtex_tissue<-gtex_tissue %>% filter(LOW_MAPABILITY<1,MAPPING_BIAS_SIM<1,GENOTYPE_WARNING<1)

gene_sym<- mapIds(org.Hs.eg.db,
                  keys=gtex_tissue$GENE_ID, #Column containing Ensembl gene ids
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multiVals="first")
gene_sym<-gene_sym[!is.na(names(gene_sym))]
gene_sym<-data.frame(symbol=unlist(gene_sym), GENE_ID=names(unlist(gene_sym)))

gtex_tissue$symbol<-gene_sym$symbol[match(gtex_tissue$GENE_ID,gene_sym$GENE_ID)]
#___________________
total_samp=length(unique(gtex_tissue$SAMPLE_ID)) #581

xx<-gtex_tissue %>% group_by(symbol,SAMPLE_ID) %>%
  summarize(total_snp=n(),sig_snp=sum(BINOM_P_ADJUSTED<0.05))%>% 
  arrange(desc(sig_snp)) %>% mutate(perc_snp=sig_snp/total_snp)

all_ase_plot<-xx %>% group_by(symbol) %>% 
  summarize(total_snp=median(total_snp),perc_snp=mean(perc_snp),n_people=sum(sig_snp>0),perc=(sum(sig_snp>0)/total_samp)*100)%>% 
  arrange(desc(perc)) 


colnames(all_ase_plot)[1]<-"gene"
gnomad<-gnomad %>%filter(canonical==TRUE) %>%  dplyr::select(gene_id,gene,lof.oe_ci.upper,transcript_type,num_coding_exons)
gnomad<-gnomad[!is.na(gnomad$transcript_type),]
all_ase_plot<-left_join(all_ase_plot,gnomad,by="gene")


xx_plot<-all_ase_plot %>%  mutate(LOEUF=cut_number(lof.oe_ci.upper,5))
quantile(all_ase_plot$perc, na.rm=T)

pdf(file="~/plot/ASE/sra_test2.pdf", width = 10, height = 6)
ggplot(all_ase_plot, aes(x = total_snp, y =lof.oe_ci.upper)) +geom_point(alpha=0.5)+ xlim(c(0,10))
ggplot(all_ase_plot, aes(x = perc_snp, y =lof.oe_ci.upper)) +geom_point(alpha=0.5)

ggplot(xx_plot, aes(x = perc_snp, y =LOEUF,fill=LOEUF)) + geom_density_ridges2()

ggplot(xx_plot, aes(x = lof.oe_ci.upper, y =n_people,fill=n_people)) + geom_density_ridges2()

dev.off()


#---------------------------------
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
exon_gr<-exons(TxDb.Hsapiens.UCSC.hg38.knownGene)


gnomad<-gnomad %>%filter(canonical==TRUE) %>%  dplyr::select(gene_id,gene,lof.oe_ci.upper,transcript_type,num_coding_exons)
gnomad<-gnomad[!is.na(gnomad$transcript_type),]
colnames(gnomad)[2]<-"symbol"

test<-gtex_tissue %>% 
  group_by(CHR,POS ,symbol,SAMPLE_ID,VARIANT_ANNOTATION) %>%
  summarize(sig_snp=sum(BINOM_P_ADJUSTED<0.05,na.rm=T)) %>% 
  ungroup()


summary_ase<-test %>% group_by(CHR,POS ,symbol,VARIANT_ANNOTATION) %>%
  summarize(n_ase=sum(sig_snp>0,na.rm=T),n_non_ase= sum(sig_snp==0,na.rm=T))
summary_ase<-left_join(summary_ase,gnomad,by="symbol")



summary_ase_gr<-makeGRangesFromDataFrame(summary_ase,seqnames="CHR",start.field ="POS",end.field = "POS")
ov<-findOverlaps(summary_ase_gr,exon_gr)
summary_ase$IsExon<-NA
summary_ase$IsExon[queryHits(ov)]<-TRUE
xx<-summary_ase %>% filter(IsExon==TRUE)

xx2<-xx %>% group_by(symbol) %>% summarize(total_ase=sum(n_ase), total_non_ase=sum(n_non_ase))
xx2<-left_join(xx2,gnomad,by="symbol")


xx2$perc_ase<-xx2$total_ase/(xx2$total_non_ase+xx2$total_ase)
xx2$perc_non_ase<-xx2$total_non_ase/(xx2$total_non_ase+xx2$total_ase)
xx2$total_snp<-(xx2$total_ase+xx2$total_non_ase)



never<-xx2 %>% filter(perc_non_ase>0.99)
dim(never)


xx2<-xx2 %>% filter(total_snp>10) %>%
  mutate(sig_ase_cut=cut_number(perc_ase,4),
         total_snp_cut=cut_number(total_snp,4),
         Loeuf_cut=cut_number(lof.oe_ci.upper,5))

pdf(file="~/plot/ASE/sra_test3.pdf", width = 8, height = 4)
ggplot(never, aes(lof.oe_ci.upper)) +
  geom_histogram()+labs(title="LOEUF for SNPs that never have sig ase")
ggplot(xx2, aes(y=sig_ase_cut,x=lof.oe_ci.upper,fill=sig_ase_cut)) +
  geom_density_ridges2()+
  geom_vline(xintercept=0.6)+
  labs(title="percent sig ase snps vs loeuf")
ggplot(xx2, aes(x=perc_ase,y=Loeuf_cut)) +
  geom_density_ridges2()+xlim(c(0,0.25))
dev.off()

xx2 %>% filter(perc_ase>.8,total_snp>30)
summary_ase %>% filter(n_ase>200)
xx2_2<-summary_ase %>% group_by(VARIANT_ANNOTATION) %>% summarize(total_ase=sum(n_ase), total_non_ase=sum(n_non_ase))
xx2<-left_join(xx2,gnomad,by="symbol")

