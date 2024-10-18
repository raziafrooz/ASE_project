library(AnnotationHub)
ah <- AnnotationHub()
ah <- query(ah, c("v26","GENCODE","Homo sapiens","GRch38")) #v26 was used in GTExV8
TxDb<-ah[["AH75155"]]


tx_38<-exons(TxDb,columns=c("TXNAME","GENEID","EXONID","CDSID"))

healthy<-final_tissue[final_tissue$dd=="No",]
i=1
sig_all<-c()

for(i in 11987:nrow(healthy)){
  print(i)
  exp_id<-healthy$experiment_acc[i]
  ase_df<-fread(metadata$ASE_path[metadata$experiment_acc == exp_id])%>% 
    mutate(alt_count=coverage-ref_count,
           log2aFC= abs((log2(alt_count+1)/ log2(ref_count+1))))

  

  #Make Granges:
  #-------------------------------
  ase_df_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="start",end.field = "start")
  
  
  ov<-findOverlaps(tx_38,ase_df_gr)
  
  ase_df$GENEID<-NA
  ase_df$GENEID[subjectHits(ov)]<-unlist(tx_38$GENEID[queryHits(ov)])
  
  ase_df$EXONID<-NA
  ase_df$EXONID[subjectHits(ov)]<-unlist(tx_38$EXONID[queryHits(ov)])
  
  
  
  # ase_df<-ase_df %>% rowwise() %>% 
  #   mutate(hap_A=max(ref_count,alt_count), hap_B=min(ref_count,alt_count))
  # 
  # 
  # exon_level<-ase_df %>% group_by(chr,GENEID,EXONID) %>%mutate(n=n()) %>% 
  #   summarize(hap_A=sum(hap_A), hap_B=sum(hap_B),
  #             ref_ratio=hap_A/(hap_A+hap_B),
  #             log2aFC= log2((hap_A+1)/(hap_B+1)))
  # 
  exon_level<-ase_df %>% filter(q_val<0.05) %>% dplyr::select(chr,start,p_val,q_val,log2aFC,GENEID,EXONID) 
  
  # exon_level[1:6,13:16]
  # exon_level[exon_level$n>530,13:16]
  # is.median<-median(exon_level$ref_ratio)
  # exon_level$p_val = apply(exon_level[,c("hap_A","hap_B")], 1, function(x) {
  #   binom.test(round(x[1],1),round((x[1]+x[2]),1),p=is.median)$p.value})
  # 
  # # perform multiple testing correction with FDR
  # exon_level$q_val = p.adjust(exon_level$p_val, method = "fdr")
  # 
  # 
  
  exon_level$sample<- exp_id
  
  
  sig_all<-rbind(sig_all,exon_level)
  
}
#fwrite(sig_all, "~/test/healthy.csv.gz")
#gene_wise<-sig_all
#sig_all<-read.csv("~/test/healthy.csv")
#  exon_wise<-sig_all
#sig_all$X<-NULL
gnomad<-read.delim("hansen_lab/dwl_files/ASE_filter/gnomad.v4.1.constraint_metrics.tsv")
gene_wise<-sig_all[-which(duplicated(sig_all$GENEID) & duplicated(sig_all$sample)),] 
x<-gene_wise  %>% group_by(GENEID) %>% summarize(n=n()) %>% arrange(desc(n))

x$n_percent<-(x$n/length(unique(sig_all$sample)))*100
gnomad<-gnomad[gnomad$canonical=="true",]
x$GENEID<-gsub("\\.\\d+", "", x$GENEID)

x$loeuf<-gnomad$lof.oe_ci.upper[match(x$GENEID,gnomad$gene_id)]

x$gene_cat[which(x$n_percent<=1)]<-"few"
x$gene_cat[which(x$n_percent>1)]<-"often"


library(ggridges)

pdf(file="~/plot/ASE/test2.pdf", width = 10, height = 6)


ggplot(x, aes(x=loeuf,y=gene_cat))+geom_density_ridges()
dev.off()
  
  
  
  