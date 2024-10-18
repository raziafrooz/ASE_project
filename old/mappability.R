library(rtracklayer)
library(GenomicRanges)
#I found a bed file including the problematic sites from ENCODE
#https://github.com/Boyle-Lab/Blacklist/tree/master
black_list<-import("~/ASE/data/hg38-blacklist.v2.bed.gz")  

#https://bismap.hoffmanlab.org
bed_graph<-import("~/k50.Unique.Mappability.bb")


sra_geno<-readRDS("~/plot/ASE/sra_qc.rds")
plot_df2[200,]
k=200
for(k in 1:length(plot_df2$sample_id)){
  print(k)
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sra_geno$sample_id[k], "_ase.rda")))
  {
    load(paste0("~/hansen_lab/ASE/test_ASE/", plot_df2$sample_id[k], "_ase.rda") ) #named ase_all or ase_df
    
    if(exists("ase_all")){
      df[k,]<-median(ase_all$ref_ratio)
      rm(ase_all)
    }else{
      df[k,]<-median(ase_df$ref_ratio)
      rm(ase_df)
    }
  }}
ase_all<-ase_df
ase_all_gr<-makeGRangesFromDataFrame(ase_all, seqnames.field = "chr", start.field = "pos", end.field = "pos")
ov<-findOverlaps(bed_graph,ase_all_gr)

x<-ase_all[unique(subjectHits(ov)),]
median(x$ref_ratio)
dim(x)
quantile(ase_all$total)

ase_all %>% filter(gene_id== "ENSG00000176155.18" ) %>% summarize(median(ref_ratio))

ov<-findOverlaps(tx_38,ase_all_gr)
quantile(ase_all$total)
ase_all$gene_id<-NA
ase_all$gene_id[subjectHits(ov)]<-as.character(tx_38$GENEID[queryHits(ov)])
length(unique(ase_all$gene_id))
ase_all %>% group_by(gene_id) %>% summarize(n=n()) %>% arrange(desc(n))
pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
ggplot(ase_all)+
  geom_point(aes(x=log2(ref), y=log2(alt)))+
  geom_point(data=ase_all[-unique(subjectHits(ov)),], aes(x=log2(ref), y=log2(alt)), color="red")
dev.off()

plot_df2[1,10:15]
ase_all[1,]
