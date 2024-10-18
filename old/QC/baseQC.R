setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(recount3)


human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  project == "SRP187978" & project_type == "data_sources"
)
rse <- create_rse(proj_info)
## Create a RangedSummarizedExperiment (RSE) object at the gene level
#rse_gene_SRP009615 <- create_rse(proj_info)
#colData(rse_gene_SRP009615)[1,30:40]
#quantile(rse_gene_SRP009615$sra.run_total_bases)
url<-locate_url(
  "SRP187978",
  "data_sources/sra",
  type = "metadata")

x <-utils::read.delim(file_retrieve(url[4], verbose = FALSE))

xx<-which(colnames(rse)%in%c("SRR8701271","SRR8701272","SRR8701273", "SRR8701274","SRR8701275","SRR8701276",
                             "SRR8701277","SRR8701278",
                             "SRR8701279","SRR8701280",
                             "SRR8701281","SRR8701282",
                             "SRR8701283","SRR8701284",
                             "SRR8701285","SRR8701286",
                             "SRR8702455","SRR8702456",
                             "SRR8702457","SRR8702458",
                             "SRR8702459","SRR8702460",
                             "SRR8702461","SRR8702462"))


tx_number<-apply(assays(rse)$raw_counts[,xx], 2, function(c)sum(c==0))

    
xx<-subset(human_projects, file_source == "gtex" & project_type == "data_sources" & project=="LIVER")
gtex <- create_rse(xx)
apply(assays(gtex)$raw_counts[,1:4], 2, function(c)sum(c!=0))


plot2<-data.frame(id="gtex",
                  indv=colData(gtex)$gtex.subjid[1:3],
                  numberOfBases= colData(gtex)[1:3,190],
                  distinct_quality_value= colData(gtex)[1:3,189],
                  tx_number= apply(assays(gtex)$raw_counts[,1:3], 2, function(c)sum(c==0)))


plot=data.frame(id=c("SRR8701271","SRR8701272","SRR8701273",
                "SRR8701274","SRR8701275","SRR8701276"), indv="indv1")

plot<-rbind(plot,data.frame(id=c("SRR8701277","SRR8701278",
                                    "SRR8701279","SRR8701280",
                                    "SRR8701281","SRR8701282",
                                    "SRR8701283","SRR8701284",
                                    "SRR8701285","SRR8701286"), indv="indv2"))

#Penn cohort1:
plot<-rbind(plot,data.frame(id= c("SRR8702455","SRR8702456",
                                    "SRR8702457","SRR8702458",
                                    "SRR8702459","SRR8702460",
                                    "SRR8702461","SRR8702462"), indv="indv3"))


plot$numberOfBases<-x$X.bases[match(plot$id,x$external_id)]
plot$distinct_quality_value<-x$X.distinct_quality_values[match(plot$id,x$external_id)]
plot$tx_number<-tx_number

pdf(file="~/plot/ASE/txQC.pdf", width = 10, height = 4)
ggplot(plot)+
  geom_jitter(aes(x=indv,y=tx_number),alpha=0.4)+
  geom_boxplot(aes(x=indv,y=tx_number),outlier.shape = NA)+
  geom_point(data=plot2,aes(x=indv,y=tx_number),alpha=0.4)+
  labs(title="number of transcripts with 0 read mapped",
       subtitle="3 GTEx liver samples vs 3 SRA indiv (with multiple runs)")
dev.off()

pdf(file="~/plot/ASE/baseQC.pdf", width = 10, height = 4)
ggplot(plot)+
  geom_point(aes(x=indv,y=numberOfBases),alpha=0.4)+
  geom_boxplot(aes(x=indv,y=numberOfBases))+
  annotate("text", label = c("#run=","6","10","8"),
    x = c(3,4,5,6), y = 100000000, size = 3, colour = "red")+
  geom_point(data=plot2,aes(x=indv,y=numberOfBases),alpha=0.4)+
  labs(title= "comparison of sequencing quality: # base count",
       subtitle= "3 samples GTEx liver vs SRP187978 project: each dot represents 1 run for indv (GTEx has 1 run)")

ggplot(plot)+
  geom_point(aes(x=indv,y=log10(numberOfBases)),alpha=0.4)+
  geom_point(data=plot2,aes(x=indv,y=log10(numberOfBases)),alpha=0.4)+
  labs(title= "comparison of sequencing quality: log10(# base count)",
       subtitle= "3 samples GTEx liver vs SRP187978 project: each dot represents 1 run for indv (GTEx has 1 run)")


ggplot(plot)+
  geom_point(aes(x=indv,y=distinct_quality_value),alpha=0.4)+
  geom_point(data=plot2,aes(x=indv,y=distinct_quality_value),alpha=0.4)+
  annotate("text", label = c("#run=","6","10","8"),
           x = c(3,4,5,6), y = 20, size = 3, colour = "red")+
  labs(title= "comparison of sequencing quality: distinct_quality_value",
       subtitle= "3 samples GTEx liver vs SRP187978 project: each dot represents 1 run for indv (GTEx has 1 run)")

dev.off()


#----------------------------------------
#Part2
#----------------------------------------

sra<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/metadata/all_SRA.csv")
sra_geno<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")

id<-which(sra_geno$sample_id %in% c("SRR8701271","SRR8701272","SRR8701273",
                                    "SRR8701274","SRR8701275","SRR8701276"))

id<-which(sra_geno$sample_id %in% c("SRR8701277","SRR8701278",
                                    "SRR8701279","SRR8701280",
                                    "SRR8701281","SRR8701282",
                                    "SRR8701283","SRR8701284",
                                    "SRR8701285","SRR8701286"))

#Penn cohort1:
id<-which(sra_geno$sample_id %in% c("SRR8702455","SRR8702456",
                                    "SRR8702457","SRR8702458",
                                    "SRR8702459","SRR8702460",
                                    "SRR8702461","SRR8702462"))
for(k in id){
  print(k)
  sample_id<-sra_geno$sample_id[k]
  study<-sra_geno$study[k]
  
  geno<-as_tibble(read.csv(sra_geno$genotypedSamples[k]))
  geno <-geno %>% filter(coverage>=8)
  
  geno$sample_id<-sample_id
  geno$study<-study
  if(k==id[1]){
    geno_all<-geno
  } else{
    geno_all<-rbind(geno_all,geno)
  }}
    
#plot<-geno_all %>% 
#  group_by(sample_id,pred_genotype) %>% 
#  summarise(number=n())
#plot$indv<- "indv1"

plot2<-geno_all %>% 
  group_by(sample_id,pred_genotype) %>% 
  summarise(number=n())
plot2$indv<- "indv3"

plot<-rbind(plot,plot2)

#---------
#Gtex
#-------
geno_met<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")
geno_met<-geno_met[geno_met$study=="Liver",]
gtex_geno<-as_tibble(readRDS(geno_met$allGenotypesOutput))
gtex_geno<-gtex_geno %>% filter(coverage>=8)

gtex_geno$sample_id<-unlist(lapply(strsplit(gtex_geno$sample_id_rep,"-"), function(xx) paste0(xx[1], "-",xx[2])))

gtex_plot<-gtex_geno %>%  filter(sample_id_rep %in% unique(gtex_geno$sample_id_rep)[1:3]) %>% 
  group_by(sample_id_rep,pred_genotype) %>% 
  summarise(number=n())
gtex_plot$indv<-gtex_plot$sample_id_rep
gtex_plot$pred_genotype<-as.integer(gtex_plot$pred_genotype)
plotXX<-rbind(plot,gtex_plot)

pdf(file="~/plot/ASE/genotype_QC.pdf", width = 10, height = 4)

ggplot(plotXX %>% filter(pred_genotype == 1))+
  geom_boxplot(aes(x=indv, y=number,color=as.factor(pred_genotype)),outlier.shape = NA)+
  geom_jitter(aes(x=indv, y=number,color=as.factor(pred_genotype)),width = 0.2, alpha=0.4)+
  labs(title="Number of SNPs called Homozygous Reference",
       color="Genotype")

ggplot(plotXX %>% filter(pred_genotype != 1))+
geom_boxplot(aes(x=indv, y=number,color=as.factor(pred_genotype)),outlier.shape = NA)+
  geom_jitter(aes(x=indv, y=number,color=as.factor(pred_genotype)), width = 0.2,alpha=0.4)+
  labs(title="Number of SNPs called (2) heterozygous and (3) homo alt",
       color="Genotype")
dev.off()

plot %>% group_by(pred_genotype, indv) %>% 
  summarise(number=mean(number))
gtex_plot %>% group_by(pred_genotype, indv) %>% 
  summarise(number=mean(number))

