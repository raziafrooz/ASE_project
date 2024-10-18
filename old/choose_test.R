setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(qvalue)
tissue<-"Liver"

load(paste0("~/hansen_lab/ASE/test_ASE/", tissue, "_ase.rda") )
true_gtex<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/true_gtex_", tissue, ".rds"))
colnames(true_gtex)[1:2]<-c("chr","pos")

sam<- unique(ase_all$sample_id)[2:10]
ase_all_1 <- ase_all %>% filter(sample_id %in% sam)
true_gtex_1 <- true_gtex %>% filter(sample_id %in% sam)

# calculate AE:
ase_all_1$ae<-abs(0.5-ase_all_1$ref_ratio)
true_gtex_1$AE<-abs(0.5-true_gtex_1$REF_RATIO)

ase_all_1$ae<-abs((ase_all_1$alt + 1)/(ase_all_1$ref + 1))
true_gtex_1$AE<-abs((true_gtex_1$ALT_COUNT + 1)/(true_gtex_1$REF_COUNT + 1))

ase_all_1$aFC_recount<-log2((ase_all_1$alt + 1)/(ase_all_1$ref + 1))
true_gtex_1$aFC_gtex<-log2((true_gtex_1$ALT_COUNT + 1)/(true_gtex_1$REF_COUNT + 1))


#pick the SNPs that are the same:
ase_same<-inner_join(true_gtex_1,ase_all_1, by=c("chr", "pos", "sample_id"))

#-------------------------
#Plot
#-------------------------

pdf(file="~/plot/ASE/hist2.pdf", width = 10, height = 6)
for(i in 2:10){
  s<-ase_same %>% filter(sample_id==unique(ase_same$sample_id)[i])
  norm_dis<-data.frame(x=seq(from = 0, to = 1, by = 0.022),y=(dbinom(x=0:45,size = 45, prob = mean(s$ref_ratio))))
  p=ggplot(s,aes(ref_ratio))+
    geom_histogram()+
    geom_histogram(data=s %>% filter(p_val<0.01),aes(ref_ratio), fill="red")+
    labs(title="Recount", subtitle = "blue: dnorm*nrow(s)")+
    geom_line(data=norm_dis,aes(x,y*nrow(s)), color="blue")#+
    #geom_freqpoly(data=norm_dis, aes(norm_dis[,1]/nrow(s)),fill="blue")
  
  print(p)
  
  norm_dis<-data.frame(x=seq(from = 0, to = 1, by = 0.022),y=(dbinom(x=0:45,size = 45, prob = mean(s$REF_RATIO))))
  p=ggplot(s,aes(REF_RATIO))+
    geom_histogram()+
    geom_histogram(data=s %>% filter(BINOM_P<0.01),aes(REF_RATIO), fill="red")+
    labs(title="GTEx", subtitle = "blue: dnorm*nrow(s)")+
    geom_line(data=norm_dis,aes(x,y*nrow(s)), color="blue")#+
    #geom_freqpoly(data=norm_dis, aes(norm_dis[,1]))
  print(p)}
dev.off()

pdf(file="~/plot/ASE/hist_overlay.pdf", width = 10, height = 6)
for(i in 2:9){
  s<-ase_same %>% filter(sample_id==unique(ase_same$sample_id)[i])
  norm_dis<-as.data.frame(rbinom(n=nrow(s),size=nrow(s), prob= mean(s$ref_ratio)))
  p=ggplot()+
    geom_histogram(data=s %>% filter(ref_ratio<=0.5),aes(ref_ratio),alpha=0.5)+
    geom_histogram(data=s %>% filter(ref_ratio>=0.5),aes(abs(ref_ratio-1)), fill="lightblue",alpha=0.5)+
    labs(title=paste0("Recount sample",i," histogram of the ref_ratio"),
         subtitle="light blue is ref_ratio >= 0.5 plot over ref_ratio<=0.5 (grey)")#+
  #geom_freqpoly(data=norm_dis, aes(norm_dis[,1]/nrow(s)),fill="blue")
  
  print(p)
  
  norm_dis<-as.data.frame(rnorm(nrow(s), mean= mean(s$REF_RATIO)))
  p=ggplot()+
    geom_histogram(data=s %>% filter(REF_RATIO<=0.5),aes(REF_RATIO),alpha=0.5)+
    geom_histogram(data=s %>% filter(REF_RATIO>=0.5),aes(abs(REF_RATIO-1)), fill="lightblue",alpha=0.5)+
    #geom_histogram(data=s %>% filter(BINOM_P<0.01),aes(REF_RATIO), fill="red")+
    labs(title=paste0("GTEx one sample",i," histogram of the ref_ratio"),
         subtitle="light blue is ref_ratio >= 0.5 plot over ref_ratio<=0.5 (grey)")#+
  #geom_freqpoly(data=norm_dis, aes(norm_dis[,1]))
  print(p)}
dev.off()
library("VGAM")

hist(pnorm(rnorm(100,0.5), 0.5))
?rnorm()
range(rnorm(100, 0.5))
# Creating the Sequence

hist(sample(gfg, 1000,replace=T))
# Plotting the beta density
plot(gfg, dbeta(gfg, 5,5), xlab="X",
     ylab = "Beta Density", type = "l",
     col = "Red")
line(gfg, dnorm(gfg, 0.5))
ggplot(as.data.frame(gfg))+
  #geom_line(aes(x=gfg,y=dbeta(gfg, 2,2)), color="red")+
  geom_point(aes(dbinom(gfg,100, 0.5)), color="blue")
gfg<-sample(gfg,100, replace = T)
dbinom(x=0:1,10, 0.5)
?dbinom

plot(dbetabinom(x=0:45,size = 45, prob =0.5,rho = 0.2))

dbinom(x=1:10,size = 80, prob = 0.5)/80
gfg = seq(0, 1, by = 0.1)
plot(dnorm(gfg,0.5)*100)
#BiocManager::install("VGAM")

x<-seq(from = 0, to = 1, by = 0.1)
dbinom(x,size = 100, prob = 0.5)
dbetabinom(x, size, prob, rho = 0, log = FALSE)


pdf(file="~/plot/ASE/beta_binomial.pdf", width = 10, height = 6)
for(i in 2:9){
  
  s<-ase_same %>% filter(sample_id==unique(ase_same$sample_id)[i])
  #norm_dis<-data.frame(x=seq(from = 0, to = 1, by = 0.022),y=(dbinom(x=0:45,size = 45, prob = mean(s$ref_ratio))))
  norm_dis<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$ref_ratio),rho = 0)))
  norm_dis.1<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$ref_ratio),rho = 0.01)))
  norm_dis.2<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$ref_ratio),rho = 0.02)))
  norm_dis.3<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$ref_ratio),rho = 0.03)))
  norm_dis.4<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$ref_ratio),rho = 0.04)))
  norm_dis.5<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$ref_ratio),rho = 0.05)))
  
  p=ggplot(s,aes(ref_ratio))+
    geom_histogram()+
    geom_line(data=norm_dis,aes(x,y*nrow(s)), color="blue")+
    geom_line(data=norm_dis.1,aes(x,y*nrow(s)), color="#2081f7")+
    geom_line(data=norm_dis.2,aes(x,y*nrow(s)), color="#18a854")+
    geom_line(data=norm_dis.3,aes(x,y*nrow(s)), color="#b38730")+
    geom_line(data=norm_dis.4,aes(x,y*nrow(s)), color="#b51105")+
    geom_line(data=norm_dis.5,aes(x,y*nrow(s)), color="salmon")+
    annotate("text", x = 0.9 , y = c(2000, 1800,1600,1400,1200,1000), label = c("rho = 0","rho = 0.01","rho = 0.02","rho = 0.03","rho = 0.04","rho = 0.05"),
    color=c("blue","#2081f7", "#18a854", "#b38730","#b51105", "salmon" ))+
    labs(title=paste0("Find best fit: Recount sample #",i," histogram of the ref_ratio"),
         subtitle="lines are betabinomial distributions with different rho values")#+
  
  print(p)
}
dev.off()



pdf(file="~/plot/ASE/beta_binomial_gtex.pdf", width = 10, height = 6)
for(i in 2:9){
  
  s<-ase_same %>% filter(sample_id==unique(ase_same$sample_id)[i])
  #norm_dis<-data.frame(x=seq(from = 0, to = 1, by = 0.022),y=(dbinom(x=0:45,size = 45, prob = mean(s$ref_ratio))))
  norm_dis<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$REF_RATIO),rho = 0)))
  norm_dis.1<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$REF_RATIO),rho = 0.01)))
  norm_dis.2<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$REF_RATIO),rho = 0.02)))
  norm_dis.3<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$REF_RATIO),rho = 0.03)))
  norm_dis.4<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$REF_RATIO),rho = 0.04)))
  norm_dis.5<-data.frame(x=seq(from = 0, to = 1, by = 0.026),y=(dbetabinom(x=0:38,size = 38, prob = mean(s$REF_RATIO),rho = 0.05)))
  
  p=ggplot(s,aes(REF_RATIO))+
    geom_histogram()+
    geom_line(data=norm_dis,aes(x,y*nrow(s)), color="blue")+
    geom_line(data=norm_dis.1,aes(x,y*nrow(s)), color="#2081f7")+
    geom_line(data=norm_dis.2,aes(x,y*nrow(s)), color="#18a854")+
    geom_line(data=norm_dis.3,aes(x,y*nrow(s)), color="#b38730")+
    geom_line(data=norm_dis.4,aes(x,y*nrow(s)), color="#b51105")+
    geom_line(data=norm_dis.5,aes(x,y*nrow(s)), color="salmon")+
    annotate("text", x = 0.9 , y = c(2000, 1800,1600,1400,1200,1000), label = c("rho = 0","rho = 0.01","rho = 0.02","rho = 0.03","rho = 0.04","rho = 0.05"),
             color=c("blue","#2081f7", "#18a854", "#b38730","#b51105", "salmon" ))+
    labs(title=paste0("Find best fit: GTEx sample #",i," histogram of the ref_ratio"),
         subtitle="lines are betabinomial distributions with different rho values")#+
  
  print(p)
}
dev.off()