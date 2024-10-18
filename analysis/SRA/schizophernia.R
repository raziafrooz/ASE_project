library(clusterProfiler)
library("org.Hs.eg.db")
library(tidyverse)
library(gwascat)
library(ggridges)
gwcat = get_cached_gwascat()

#xx<-unique(gwcat$`DISEASE/TRAIT`)
#<-gwcat[which(str_detect(xx, "schizophrenia")==TRUE),]
xx<-gwcat %>% filter(str_detect(`DISEASE/TRAIT`,"schizophrenia"))
#unique(xx$`REPORTED GENE(S)`)

sig_all<-read.csv("~/test/schizophernia.csv")
sig_all$X<-NULL
sig_all<-sig_all %>% mutate(log2aFC= (log2(ALT_COUNT+1)/log2(REF_COUNT+1))) %>% unique()
# Get LOEUF scores:
gnomad<-read.delim("hansen_lab/dwl_files/ASE_filter/gnomad.v4.1.constraint_metrics.tsv")
gnomad<-gnomad[gnomad$canonical=="true",]
#length(unique(sig_all$sample_id))

#xx_plot<-sig_all[which(sig_all$adj.pval<0.05),] %>% group_by(chr,start,GENE_ID) %>% summarize(n=length(unique(sample_id)))%>%  arrange(desc(n))
xx_plot<-sig_all[which(sig_all$adj.pval<0.05),] %>% 
  group_by(GENE_ID) %>% 
  summarize(n=length(unique(sample_id)),m_aFC=median(abs(log2aFC))) %>%
  arrange(desc(n))

#ENSG00000111371,ENSG00000126777
# xx_plot<-sig_all[which(sig_all$adj.pval<0.05),] %>% group_by(chr,GENE_ID,sample_id) %>% summarize(n=n()) %>%  arrange(desc(n))
# xx_plot<-xx_plot %>% group_by(GENE_ID,chr) %>% summarize(n=n()) %>%  arrange(desc(n))
# xx_plot<-xx_plot[xx_plot$n>=10,]

xx_plot$loeuf<-gnomad$lof.oe_ci.upper[match(xx_plot$GENE_ID,gnomad$gene_id)]
# xx_plot$loeuf_ntile<-as.factor(ntile(xx_plot$loeuf,10))
# xx_plot$n_ntile<-as.factor(ntile(xx_plot$n,10))

xx_plot$symbol = mapIds(org.Hs.eg.db,
                        keys=xx_plot$GENE_ID, #Column containing Ensembl gene ids
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

all_n<-c()
for(gene_sym in unique(xx_plot$symbol) ){
  print(gene_sym)
symb_n<-gwcat %>% filter(str_detect(MAPPED_GENE,gene_sym)) %>% summarize(n=length(unique(LINK)))
symb_n$symbol<-gene_sym

all_n<-rbind(all_n,symb_n)
}

colnames(all_n)[1]<-"GWAS_n"
xx_plot<-left_join(xx_plot,all_n)

xx_plot[1:7,]

top_genes<-xx_plot[which(xx_plot$n>50),]
#--------------------
xx_plot<-x[!is.na(x$GENEID),]
xx_plot$GENE_ID<-x$GENEID[!is.na(x$GENEID)]


xx_plot$entrez = mapIds(org.Hs.eg.db,
                     keys=xx_plot$GENE_ID, #Column containing Ensembl gene ids
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

# 
# library(DOSE)
# 
# x <- enrichDO(gene          = xx_plot$entrez,
#               ont           = "DO",
#               pvalueCutoff  = 0.05,
#               pAdjustMethod = "BH",
#               universe      = xx_plot$chr,
#               minGSSize     = 5,
#               maxGSSize     = 500,
#               qvalueCutoff  = 0.05,
#               readable      = FALSE)

# y <- gseDO(xx_plot$entrez,
#            minGSSize     = 120,
#            pvalueCutoff  = 0.2,
#            pAdjustMethod = "BH",
#            verbose       = FALSE)

dgn <- enrichDGN(xx_plot$entrez)

#xx_plot$entrez<-do.call(rbind.data.frame,xx_plot$entrez)


# xx_plot_genes <- bitr(xx_plot$entrez, fromType = "ENTREZID",
#                       toType = c("ENSEMBL", "SYMBOL"),
#                       OrgDb = org.Hs.eg.db)


ggo <- groupGO(gene     = xx_plot$entrez,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

ego <- enrichGO(gene          = xx_plot$entrez,
                universe      = xx_plot$chr,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1, minGSSize = 3, maxGSSize = 3000,
                readable      = TRUE, pool = TRUE)




head(ego, 3)[1:2,]


pdf(file="~/plot/ASE/schizophernia.pdf", width = 7.5, height = 9)

dotplot(ego, showCategory=30) + ggtitle("schizophernia Gene enrichment analysis")
dotplot(dgn, showCategory=30) + ggtitle("schizophernia Over-representation analysis for the disease gene network")

# barplot(ego, showCategory = 20) 
# mutate(ego, qscore = -log(p.adjust, base = 10)) %>% 
#   barplot(x = "qscore")

dev.off()



#--------------------
xx_plot$ASE_pop_percent<-as.factor(ntile(xx_plot$n/length(unique(sig_all$sample_id)),3))

xx_plot[1:18,]
pdf(file="~/plot/ASE/test.pdf", width = 10, height = 6)


ggplot(xx_plot, aes(x=loeuf,y=ASE_pop_percent, fill=ASE_pop_percent))+
  geom_density_ridges(alpha=0.4)+
  labs(title=paste0("loeuf in ",nrow(xx_plot), " healthy SRA indv"))+
  geom_vline(xintercept = 0.6, color="red")
dev.off()


