setwd("~/ASE/")
library(tidyverse)
library(ggplot2)
library(ggridges)
library(data.table)


recount3_metadata<-read_tsv("/dcs04/hansen/data/recount_genotype/PCA/SRA/Recount3_metadata.tsv")
battle<-read_tsv("data/battle.tsv")

# colnames(recount3_metadata)[1:5]
# recount3_metadata[10:20,54:58]
# quantile(recount3_metadata$bc_frag.count)

#things to look at:
#library_selection, library_layout, sample_bases, size,run_total_bases,
#bc_auc.all_reads_all_bases, bc_auc.unique_reads_all_bases,bc_frag.count,bc_frag.mode_length
#data$bc_auc.unique_reads_all_bases[1:20]/data$bc_auc.all_reads_all_bases[1:20]
#bases

recount3_metadata$cancer<-battle$cancer[match(recount3_metadata$experiment_acc,battle$sample)]
recount3_metadata_nonCancer<-recount3_metadata[-which(recount3_metadata$cancer == "cancer"),]

paired<-recount3_metadata_nonCancer[recount3_metadata_nonCancer$library_layout == "paired",]
single<-recount3_metadata_nonCancer[recount3_metadata_nonCancer$library_layout == "single",]


sra_subset<-readRDS("~/test/m.rds")
sra_subset$cancer<-battle$cancer[match(sra_subset$study,battle$study)]

sra_subset_nonCancer<-sra_subset[-which(sra_subset$cancer=="cancer"),]
sra_subset_Cancer<-sra_subset[which(sra_subset$cancer=="cancer"),]

#sra_subset<-readRDS("~/plot/ASE/sra_subset.rds")
#sra_subset$disease.category<-battle$disease.category[match(sra_subset$study,battle$study)]
dim(paired)

pdf(file="~/plot/ASE/test/deming.pdf", width = 10, height = 6) #library_selection_NOTcDNA_RTPCT.pdf

data<-paired[which(sra_subset_nonCancer$study %in% paired$study),]
for(k in 8:15){
  sample_id<-data$external_id[k]
  study<-data$study[k]
  cc<-data$cancer[k]
  #dd<-data$disease.category[k]
  ll<-data$library_selection[k]
  bb<-round(data$bc_auc.unique_reads_all_bases[k]/data$bc_auc.all_reads_all_bases[k],2)
  bases<-round(data[k,'#bases']/10^8,2)
  print(k)
  
  
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
  load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
  if(exists("ase_all")){
    ase_df<-ase_all
    rm(ase_all)
  } else{
    #data$num_snp[k]<- nrow(ase_df)
    
    
    xx_lm<-round(summary(lm(log10(alt)~log10(ref), data=ase_df))$r.squared,2)
    xx_sigma<-round(summary(lm(log10(alt)~log10(ref), data=ase_df))$sigma,2)
    lm_deming_inter<-deming(log10(alt)~log10(ref), data=ase_df)$coefficients[[1]]
    lm_deming_slope<-deming(log10(alt)~log10(ref), data=ase_df)$coefficients[[2]]
    
    p0=ggplot(ase_df, aes(y=log10(alt), x=log10(ref)))+
      geom_point(alpha=0.4)+
      geom_smooth(method = lm, se = FALSE, color="blue")+
      geom_abline(intercept = lm_deming_inter, slope = lm_deming_slope,color="green")+
      geom_abline(intercept = 0, slope = 1,color="red")+
      geom_abline(intercept = 0.25, slope = 1,color="salmon")+
      geom_abline(intercept = -0.25, slope = 1,color="salmon")+
      labs(title= paste0(sample_id,"-", data$library_layout[k]),
           subtitle= paste0(study,"-",round(median(ase_df$ref_ratio),2),
                            "-nrow:", nrow(ase_df),"-",cc,"-", ll, "-unique:",bb,
                            "-#bases:", bases))+
    xlim(c(0,4))+
      ylim(c(0,4))+
      annotate("text", x = 0.5, y = 3.5,
               label = paste0("R^ 2 = ",xx_lm))+
      annotate("text", x = 0.5, y = 3,
               label = paste0("sigma = ",xx_sigma))
    print(p0)
    
  }
  
  rm(ase_df)
  
}}
dev.off()


install.packages("MethComp")
library(MethComp)
m1 <- Deming(ase_df$alt,ase_df$ref)
m1
paired<-recount3_metadata[recount3_metadata$library_layout == "paired",]
dim(paired)
lm_deming_inter<-deming(alt~ref, data=ase_df)$coefficients[[1]]
lm_deming_slope<-deming(alt~ref, data=ase_df)$coefficients[[2]]

#Filter1:
#filter based on % uniqueness
#------------------
paired$unique<-paired$bc_auc.unique_reads_all_bases/paired$bc_auc.all_reads_all_bases
paired1<-paired[which(paired$unique>=0.7),]
dim(paired1)

#Filter2:
#Filter the samples that have total base number higher than 50e8
#-------------------
paired2<-paired1[which(paired1$'#bases'<50e8),]
dim(paired2)

#Filter3:
#Filter the samples with extreme 
#-------------------
paired2<-paired1[which(paired1$'#bases'<50e8),]
dim(paired2)





#-----------------------------------------------------------------
pdf(file="~/plot/ASE/test/MA1.pdf", width = 10, height = 6) #library_selection_NOTcDNA_RTPCT.pdf

data<-paired[which(sra_subset_nonCancer$study %in% paired$study),]
for(k in 1:20){
  sample_id<-data$external_id[k]
  study<-data$study[k]
  cc<-data$cancer[k]
  #dd<-data$disease.category[k]
  ll<-data$library_selection[k]
  bb<-round(data$bc_auc.unique_reads_all_bases[k]/data$bc_auc.all_reads_all_bases[k],2)
  bases<-round(data[k,'#bases']/10^8,2)
  print(k)
  
  
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
    if(exists("ase_all")){
      ase_df<-ase_all
      rm(ase_all)
    } else{
      
ase_df<- ase_df %>% mutate(ratio=log2(alt) - log2(ref),
                               mean=(log2(alt) + log2(ref))/2)

ase_df$ntile<-ntile(ase_df$mean,20)
median(ase_df$ratio)
positive_sd<-ase_df %>% group_by(ntile) %>% mutate(ratio_mean=mean(ratio)) %>% 
  filter(ratio>ratio_mean) %>% summarize(max=max(mean),sd=sd(ratio) )
negative_sd<-ase_df %>% group_by(ntile) %>% mutate(ratio_mean=mean(ratio)) %>% 
  filter(ratio<ratio_mean) %>% summarize(max=max(mean),sd=sd(ratio))
mean<-ase_df %>% group_by(ntile) %>% summarize(max=max(mean),mean=mean(ratio) )

positive_sd$sd<-positive_sd$sd+mean$mean
negative_sd$sd<-(-(negative_sd$sd+mean$mean))

p0=ggplot(ase_df, aes(y=ratio, x=mean))+
  geom_jitter(alpha=0.4)+
  geom_line(data=positive_sd, aes(x=max, y=sd, group=1), color="red")+
  geom_line(data=negative_sd, aes(x=max, y=sd, group=1), color="red")+
  geom_line(data=mean, aes(x=max, y=mean, group=1), color="green")+
  geom_hline(yintercept = 0, color="black")+
  labs(title= paste0(sample_id,"-", data$library_layout[k]),
       subtitle= paste0(study,"-",round(median(ase_df$ref_ratio),2),
                        "-nrow:", nrow(ase_df),"-",cc,"-", ll, "-unique:",bb,
                        "-#bases:", bases))+
annotate("text", x = 3.5, y = 2.5,
         label = max(positive_sd$sd))+
  annotate("text", x = 3.5, y = 2,
           label = max(negative_sd$sd))
print(p0)

p0=ggplot(ase_df, aes(ref_ratio))+
  geom_histogram()
print(p0)
    }
    
    rm(ase_df)
    
  }}
dev.off()



library(rstatix)
sample_id<-data$external_id[17]
  load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
  
ase_df%>% sample_n(500, replace = FALSE) %>% shapiro_test(ref_ratio)
library(mecor)
library(deming)
library(MethComp)
x_lm<-lm(alt~ref, data=ase_df)
x_lm2<-Deming(y=log2(ase_df$alt),x=log2(ase_df$ref),boot=TRUE)
rc_fit <- 
  mecor(formula = alt ~ MeasError(substitute = ref, reference = alt),
        data = ase_df,
        method = "standard", # defaults to "standard"
        B = 0)

summary(rc_fit)

summary(x_lm)
x_lm2$coefficients[[2]]
x_lm$coefficients[[1]]
median()

#-----------------------------------------------------------------
data<-paired[which(sra_subset_nonCancer$study %in% paired$study),]
k=7
pdf(file="~/plot/ASE/test/MA_nonlog1.pdf", width = 10, height = 6) #library_selection_NOTcDNA_RTPCT.pdf
for(k in 7:12){
  sample_id<-data$external_id[k]
  study<-data$study[k]
  cc<-data$cancer[k]
  #dd<-data$disease.category[k]
  ll<-data$library_selection[k]
  bb<-round(data$bc_auc.unique_reads_all_bases[k]/data$bc_auc.all_reads_all_bases[k],2)
  bases<-round(data[k,'#bases']/10^8,2)
  print(k)
  
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
    if(exists("ase_all")){
      ase_df<-ase_all
      rm(ase_all)
    } else{
      
    
    ase_df<- ase_df %>% mutate(ratio=log10(alt) - log10(ref),
                               mean=(log10(alt) + log10(ref))/2)
    
  
    
    
    
    
    ase_df$mean_20tile<-ntile(ase_df$mean,30)
    
    positive<-ase_df %>% group_by(mean_20tile) %>%
      mutate(median_ntile=median(ratio)) %>% ungroup() %>% 
      filter(ratio>median_ntile) %>% group_by(mean_20tile) %>%
      summarise(enframe(quantile(ratio, 0.75), "quantile", "positive")) 
      
    negative<-ase_df %>% group_by(mean_20tile) %>%
      mutate(median_ntile=median(ratio)) %>% ungroup() %>% 
      filter(ratio<median_ntile) %>% group_by(mean_20tile) %>%
      summarise(enframe(quantile(ratio, 0.75), "quantile", "negative"))
    
    
    quantile_lines<-ase_df %>% group_by(mean_20tile) %>%
      summarize(median_ntile=median(ratio), max_mean=max(mean))%>% ungroup() %>% 
      left_join(positive) %>% left_join(negative)
    
    quantile_lines$positive_invert<-(-quantile_lines$positive)
    
    
    p0=ggplot(ase_df, aes(y=ratio, x=mean))+
      geom_point(alpha=0.4)+
      geom_hline(yintercept = 0, color="black")+
      #geom_quantile(method = "rqss", color="red",lambda = 0.01,alpha=0.5, quantiles = q10)+
      labs(title= paste0(sample_id,"-", data$library_layout[k]),
           subtitle= paste0(study,"-",round(median(ase_df$ref_ratio),2),
                            "-nrow:", nrow(ase_df),"-",cc,"-", ll, "-unique:",bb,
                            "-#bases:", bases))+
      geom_line(data=quantile_lines, aes(x=max_mean, y=positive, group=1), color="red")+
      geom_line(data=quantile_lines, aes(x=max_mean, y=negative, group=1), color="red")+
      geom_line(data=quantile_lines, aes(x=max_mean, y=positive_invert, group=1), color="purple")+
      geom_line(data=quantile_lines, aes(x=max_mean, y=median_ntile, group=1), color="green")
      #xlim(c(0,quantile(ase_df$mean,0.95)[[1]]))+
      #ylim(c(-100,100))
    
    print(p0)
    rm(ase_df)
}}}
    dev.off()
    q10 <- seq(0.05, 0.95, by = 0.05)
    quantile(ase_df$ratio)
      # 
      # 
      # 
# 
# 
# 
# 
# 
# 

    pdf(file="~/plot/ASE/test/adjust2.pdf", width = 10, height = 6) #library_selection_NOTcDNA_RTPCT.pdf
    for(k in 7:19){
      sample_id<-data$external_id[k]
      study<-data$study[k]
      cc<-data$cancer[k]
      #dd<-data$disease.category[k]
      ll<-data$library_selection[k]
      bb<-round(data$bc_auc.unique_reads_all_bases[k]/data$bc_auc.all_reads_all_bases[k],2)
      bases<-round(data[k,'#bases']/10^8,2)
      print(k)
      
      if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
        load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
        if(exists("ase_all")){
          ase_df<-ase_all
          rm(ase_all)
        } else{
          
        if(median(ase_df$ref_ratio)<0.5){
          
          adj<-round(ase_df$alt *0.1,0)
          ase_df$new_ref<-ase_df$ref+adj
          ase_df$new_alt<-ase_df$alt-adj
        }else if(median(ase_df$ref_ratio)>0.51){
          
          adj<-round(ase_df$ref *0.1,0)
          ase_df$new_ref<-ase_df$ref-adj
          ase_df$new_alt<-ase_df$alt+adj
          
        }else{
          ase_df$new_ref<-ase_df$ref
          ase_df$new_alt<-ase_df$alt
        }
          
          p0=ggplot(ase_df, aes(y=log10(new_alt), x=log10(new_ref)))+
            geom_point(alpha=0.2, color="black")+
            #geom_point(data=ase_df, aes(y=log10(alt), x=log10(ref)), alpha=0.1)+
            labs(title= paste0(sample_id,"-", data$library_layout[k]),
                 subtitle= paste0(study,"-","Old ref_ratio:",round(median(ase_df$ref_ratio),2),"-",
                                  "new_ref_ratio:", 
                                  round(median(ase_df$new_ref/ase_df$total),2)))+
            geom_abline(intercept = 0, slope = 1,color="red")
          print(p0)
          
        }}}
    dev.off()
          
          
    
    
    
    
    k=7
    
    
    pdf(file="~/plot/ASE/test/test.pdf", width = 10, height = 6) #library_selection_NOTcDNA_RTPCT.pdf
    seq_mean=seq(0,4.6,by=0.5)
    for(k in 7:12){
      sample_id<-data$external_id[k]
      study<-data$study[k]
      print(k)
      
      if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
        load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
        if(exists("ase_all")){
          ase_df<-ase_all
          rm(ase_all)
        } else{
          print("all good")
          }
        
          ase_df<- ase_df %>% mutate(ratio=log10(alt) - log10(ref),
                                     mean=(log10(alt) + log10(ref))/2)
          
          
          
          ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)
          
          q_line<-ase_df %>% group_by(mean_20tile) %>%
            summarise(median_ratio=median(ratio),
                      enframe(quantile(ratio, c(0.05,0.95)), "quantile", "ratio_q"),
                      min=min(mean),max=max(mean)) %>% 
            mutate(q=case_when(quantile== '5%' ~ "low",
                               quantile== '95%' ~ "high" ))
          
          positive_invert<- q_line[which(q_line$q=="high"),] %>% 
            mutate(dist= ratio_q - median_ratio, 
                   invert= median_ratio-dist)
            
          line_segment<-ase_df %>% group_by(mean_20tile) %>%
            summarise(min=min(mean),max=max(mean))
          
          
          p0=ggplot(ase_df, aes(y=ratio, x=mean))+
            geom_point(alpha=0.4)+
            geom_hline(yintercept = 0, color="black")+
            labs(title= paste0(sample_id,"-", data$library_layout[k]),
                 subtitle= paste0(study,"-","Old ref_ratio:",round(median(ase_df$ref_ratio),2)))
          
          p1=p0+
             geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=ratio_q, group=1), color="red")+
            geom_line(data=q_line[which(q_line$q=="low"),], aes(x=max, y=ratio_q, group=1), color="red")+
             geom_line(data=positive_invert, aes(x=max, y=invert, group=1), color="purple")+
            geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=median_ratio, group=1), color="green")+
            geom_segment(data = line_segment, mapping = aes(x=min, y=1, xend=max, yend=1), inherit.aes = FALSE, color="blue",
                         arrow = arrow(length = unit(0.1,"cm")))
          print(p1)
          p2=p0+
            geom_quantile(quantiles = c(0.05, 0.5, 0.95),
                          color="orange")
          
          print(p2)
          
          
          p0=ggplot(ase_df, aes(y=log10(alt), x=log10(ref)))+
            geom_point(alpha=0.4)+
            geom_hline(yintercept = 0, color="black")+
            labs(title= paste0(sample_id,"-", data$library_layout[k]),
                 subtitle= paste0(study,"-","Old ref_ratio:",round(median(ase_df$ref_ratio),2)))+
            geom_quantile(quantiles = c(0.05, 0.5, 0.95),
                          color="orange", )
          
          print(p0)
          
          
          rm(ase_df)
        }}
    dev.off()
    
    
    
    
    