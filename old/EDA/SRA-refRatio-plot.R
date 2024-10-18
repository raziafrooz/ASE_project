setwd("~/ASE/")
library(tidyverse)
library(ggridges)
plot1<-readRDS("~/ASE/data/use_allSRA_QC.rds")
#df_test<-readRDS("~/plot/ASE/test.rds")
#df<-readRDS("~/ASE/data/allSRA_QC.rds")

#df_nor<-readRDS("~/plot/ASE/sra_qc_normal.rds")
#plot1<-inner_join(df,df_test)
#plot1<-bind_rows(plot1, df_nor)
#x<-readRDS("~/plot/ASE/coverage_test.rds")
#plot1<-inner_join(plot1,x)
#saveRDS(plot1, file="~/ASE/data/use_allSRA_QC.rds")

#---------------------------------
#df<-readRDS("~/plot/ASE/sra_qc.rds")
#colnames(df)

plot1$dis<-"Medium"
plot1$dis[plot1$ref_ratio<0.48]<-"Low"
plot1$dis[plot1$ref_ratio>0.52]<-"High"


plot2<-plot1 %>% mutate(gene_interv= cut_number(num_genes, 5),
                        ref_int= cut_number(ref_ratio, 5),
                        num_int= cut_number(avg_len, 3),
                        geno_int= cut_number(geno_cov, 6),
                        hom_ref_int= cut_number(homo_ref_cov, 6),
                        het_int= cut_number(het_cov, 5),
                        len_int= cut_number(avg_len, 3),
                        nSNP_int=cut_number(nSNP, 5)) 
plot2<-plot2 %>% group_by(dis) %>% mutate(ref_int=cut_number(ref_ratio, 2))
plot2$ref_int<-factor(plot2$ref_int, levels= c("[0.176,0.437]","(0.437,0.48]",
                                              "[0.48,0.5]","(0.5,0.52]",
                                              "[0.52,0.529]","(0.529,0.6]"))
quantile(plot2$num_genes, na.rm=T)

pdf(file="~/plot/ASE/sra_nSNP.pdf", width = 10, height = 4)

ggplot(plot2,aes(y=ref_int, x=num_genes, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,1.5e+4))

ggplot(plot2,aes(y=ref_int, x=nSNP, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,4e+4))

ggplot(plot2,aes(y=ref_int, x=nSNP/num_genes, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,10))

ggplot(plot2,aes(y=ref_int, x=ase_sig, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,8e+3))

ggplot(plot2,aes(y=ref_int, x=(ase_sig/nSNP)*100, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,40))
dev.off()

pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
ggplot(plot2,aes(y=len_int, x=num_genes))+
  geom_density_ridges(alpha=0.5)

ggplot(plot2,aes(y=len_int, x=nSNP))+
  geom_density_ridges(alpha=0.5)

ggplot(plot2,aes(y=len_int, x=ase_sig))+
  geom_density_ridges(alpha=0.5)

ggplot(plot2,aes(y=len_int, x=ase_sig))+
  geom_density_ridges(alpha=0.5)

dev.off()
colnames(plot2)


pdf(file="~/plot/ASE/sra_coverage.pdf", width = 10, height = 4)

ggplot(plot2,aes(y=ref_int, x=geno_cov, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,2e+8))

ggplot(plot2,aes(y=ref_int, x=het_cov/nSNP, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,300))

ggplot(plot2,aes(y=ref_int, x=het_cov, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,5e+6))

ggplot(plot2,aes(y=ref_int, x=homo_ref_cov, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,2e+8))

ggplot(plot2,aes(y=ref_int, x=homo_ref_cov/het_cov, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,150))
dev.off()

pdf(file="~/plot/ASE/SRA-QC.pdf", width = 10, height = 4)
ggplot(plot2,aes(y=ref_int, x=avg_len, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,200))+
  labs(title="Average sequencing length between 3 distibution",
       subtitle="low = ref_ratio<0.48, high= ref_ratio>0.52")

ggplot(plot2,aes(y=ref_int, x=numberOfBases, fill=dis))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,1e+10))+
  labs(title="Number of bases between 3 distibution",
       subtitle="low = ref_ratio<0.48, high= ref_ratio>0.52")
dev.off()

plot1<-plot1 %>% mutate(len_int2= cut_number(avg_len, 3))
plot1[which(is.na(plot1$len_int2))[1:3],]
pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
ggplot(plot1 %>% filter(study %in% unique(plot1$study)[1:100]),aes(y=ref_ratio, x=study))+
  geom_point()

ggplot(plot,aes(y=ref_int, x=nSNP))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,1e+8))
ggplot(plot1,aes(y=dis, x=het_cov))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,1e+6))
ggplot(plot1,aes(y=dis, x=het_cov, fill=len_int2))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,1e+6))
ggplot(plot1,aes(y=dis, x=homo_ref_cov, fill=len_int2))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,1e+8))
ggplot(plot1,aes(y=dis, x=homo_alt_cov, fill=len_int2))+
  geom_density_ridges(alpha=0.5)+
  xlim(c(0,1e+7))
dev.off()

quantile(plot1$num_genes[which(plot1$ref_ratio<=0.4)])

