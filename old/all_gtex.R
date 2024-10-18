
setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

genotyped_individuals <- readLines("/dcl01/hansen/data/gtex_private/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.header.txt")

geno_met<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")
geno_met$group<-"testing"
geno_met2<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Training_metadata.csv")
colnames(geno_met2)
g<-readRDS(geno_met2$prior)
geno_met<-geno_met[,c("study","allGenotypesOutput")]


bad_snp<-read.table("~/plot/ASE/hg38_haplo_count_blacklist.chr.bed", sep="\t")
bad_snp_gr<-makeGRangesFromDataFrame(bad_snp,seqnames="V1",start.field ="V2",end.field = "V3")


geno_met[1,]
colnames(geno_met)

k=28
for(k in 1:length(geno_met$study)){
  #sample_id<-geno_met$sample_id_rep[k]
  study<-geno_met$study[k]
  print(study)
  geno_dir<-geno_met$allGenotypesOutput[k]
  if(!file.exists(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rda")))
  {
    
    ase_df<-as_tibble(readRDS(geno_dir))%>% 
      filter(pred_genotype==2, coverage>=8)
    
    if(nrow(ase_df)>10){
      ase_df <-ase_df %>%
        mutate(alt=(sqrt(2^((2*S)-M))-1),ref=((2^M)*(alt+1))-1) %>% 
        select(!c(M,S,pred_genotype))
      colnames(ase_df)[c(2,4)]<-c("pos","total")
      #ase_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, ".rds"))
      ase_df$ref_ratio<- ase_df$ref/ase_df$total
      ase_df$aFC<-log2((ase_df$alt+1)/(ase_df$ref+1))
      
      #Make Granges:
      ase_filt_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="pos",end.field = "pos")
      
      #Remove blacklist snps from our granges:
      ov<-findOverlaps(ase_filt_gr,bad_snp_gr)
      print(paste0("number of blacklist: ", length(unique(queryHits(ov)))))
      if(length(unique(queryHits(ov)))>0){
        ase_df<-ase_df[-unique(queryHits(ov)),]
        
        #obtain the median of the ref ratio to be used as the null p value
        is.median<-median(ase_df$ref_ratio)
        print(paste0("median is ",is.median))
        #ase_df$p_val = apply(ase_df[,c("ref","alt")], 1, function(x) binom.test(x[1],(x[1]+x[2]),p=is.median)$p.value)
        # perform multiple testing correction with FDR
        #ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")
        
        save(ase_df, file=paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rda") )
      }
    }}}
k=2

#------------------------------------------------------
#Get the true ASE hits from gtex 
#------------------------------------------------------
met<-read.csv("data/GTEx_metadata.csv")
tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]


i=4
for(i in 1:nrow(tissues_names)){
  print(i)
  print(tissues_names$full_name[i])
  abb<-tissues_names$abb[i]
  tissue<-tissues_names$full_name[i]
  path<-tissues_names$file_name[i]
  gtex<-"/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"
  gtex_tissue<- fread(paste0(gtex,path))
  
  load(paste0("~/hansen_lab/ASE/test_ASE/", tissue, "_ase_MS.rda"))
  ase_df$sample_id<-sapply(strsplit(ase_df$sample_id_rep,"-"), function(xx){ paste0(xx[1],"-",xx[2]) })
  gtex_tissue<-as_tibble(gtex_tissue[which(gtex_tissue$SUBJECT_ID %in% ase_df$sample_id),])
  
  for(mm in 1:length(unique(ase_df$sample_id))){
    m<-unique(ase_df$sample_id)[mm]
    t1<-  gtex_tissue %>% filter(SUBJECT_ID==m)
    t1.1 <-ase_df %>% filter(sample_id==m)
    
  gr_ase<-makeGRangesFromDataFrame(t1.1,
                           keep.extra.columns=FALSE,
                           seqnames.field="chr",
                           start.field="pos",
                           end.field="pos")
  gr_gtex<-makeGRangesFromDataFrame(t1,
                                   keep.extra.columns=FALSE,
                                   seqnames.field="CHR",
                                   start.field="POS",
                                   end.field="POS")
  ov<-findOverlaps(gr_ase,gr_gtex)
  
t1 <- t1[subjectHits(ov),]

if(mm==1){
  true_gtex<-t1
  } else {
    true_gtex<-rbind(true_gtex,t1)
  }
}
saveRDS(true_gtex,file = paste0("~/hansen_lab/ASE/test_ASE/true_gtex_", tissue, ".rds"))
}

dd<-true_gtex %>% group_by(SUBJECT_ID) %>% summarize(ref_ratio=median(REF_RATIO))
df$sample_id<-sapply(strsplit(df $sample_id_rep,"-"), function(xx){ paste0(xx[1],"-",xx[2]) })
df$true_ref_ratio<-dd$ref_ratio[match(df$sample_id,dd$SUBJECT_ID)]
pdf(file="~/plot/ASE/gtex_true_ours.pdf", width = 10, height = 4)
ggplot(df)+
  geom_point(aes(x=true_ref_ratio, y=ref_ratio),alpha=0.4)
dev.off()

#-------------------------------
#Plot
#--------------------------------
for(k in 1:length(geno_met$study)){
  study<-geno_met$study[k]
  print(study)
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rda")))
  {
    load(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rda"))
    
    df<-ase_df %>% group_by(sample_id_rep) %>% summarize(ref_ratio=median(ref_ratio))
    df$tissue<-study
    if(k==1){
      ase_gtex<-df
    }
    ase_gtex<-rbind(ase_gtex,df)
  }
}
#saveRDS(ase_gtex, file="data/gtex_plot.rds")
rr=0.53
g<-ase_gtex[which(ase_gtex$ref_ratio>=rr),]
g$sample_id<-unlist(lapply(strsplit(g$sample_id_rep, "-"), function(xx) {paste0(xx[[1]], "-", xx[[2]])}))
unique(g$sample_id)
pdf(file="~/plot/ASE/gtex_all.pdf", width = 10, height = 4)
ggplot(data=ase_gtex[which(ase_gtex$ref_ratio>rr),], aes(x=tissue,y=ref_ratio))+
  geom_point(alpha=0.5)+
  geom_point(data= g, 
             aes(x=tissue,y=ref_ratio, color=sample_id),alpha=0.5)+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  labs(title="Ref_ratio across 54 tissues in GTEx (test set = 4076 samples)",
       subtitle = "Each dot is the median of ref_ratio in one sample")

ggplot(data=ase_gtex[which(ase_gtex$ref_ratio>rr),], aes(x=tissue,y=ref_ratio))+
  geom_point(alpha=0.5)+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  labs(title="Ref_ratio across 54 tissues in GTEx (without outliers)",
       subtitle = "Each dot is the median of ref_ratio in one sample")
dev.off()
dim(ase_gtex)

