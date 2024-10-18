#geuvadis ERP001942 
library(quantreg)

geu<-recount3_metadata[which(recount3_metadata$study=="ERP001942"),]
bad_samples<-c("ERR188390","ERR204881","ERR204916","ERR204940",
               "ERR204843","ERR204897","ERR205016","ERR188328")
# 
# pdf(file="~/plot/ASE/test/geu.pdf", width = 10, height = 6) 
# for(k in 1:nrow(geu)){
#   sample_id<-geu$external_id[k]
#   study<-geu$study[k]
#   print(k)
#   
#   if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
#     load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
#     if(exists("ase_all")){
#       ase_df<-ase_all
#       rm(ase_all)
#     } else{
#       print("ase found")}
#     
#     p0=ggplot(ase_df, aes(y=log10(alt), x=log10(ref)))+
#       geom_point(alpha=0.4)+
#       labs(title= paste0(sample_id,"-", data$library_layout[k]),
#            subtitle= paste0(study,"-","Old ref_ratio:",round(median(ase_df$ref_ratio),2)))
#     print(p0)
#     
#     
#     rm(ase_df)
#   }}
# dev.off()   
# 


# 
# pdf(file="~/plot/ASE/test/quantile_adjust_geu.pdf", width = 10, height = 6) 
# for(k in 1:50){
#   sample_id<-geu$external_id[k]
#    study<-geu$study[k]
#   print(k)
#   
#   if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
#     load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
#     if(exists("ase_all")){
#       ase_df<-ase_all
#       rm(ase_all)
#     } else{
#       print("ase found")}
#     ase_df<- ase_df %>% mutate(ratio=log10(alt) - log10(ref),
#                                mean=(log10(alt) + log10(ref))/2)
#     
#    
#     ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)
#     
#     q_line<-ase_df %>% group_by(mean_20tile) %>%
#       summarise(median_ratio=median(ratio),
#                 enframe(quantile(ratio, c(0.05,0.95)), "quantile", "ratio_q"),
#                 min=min(mean),max=max(mean)) %>% 
#       mutate(q=case_when(quantile== '5%' ~ "low",
#                          quantile== '95%' ~ "high" ))
#     
#     positive_invert<- q_line[which(q_line$q=="high"),] %>% 
#       mutate(dist= ratio_q - median_ratio, 
#              invert= median_ratio-dist)
#     
#     line_segment<-ase_df %>% group_by(mean_20tile) %>%
#       summarise(min=min(mean),max=max(mean))
#     
#     
#     p0=ggplot(ase_df, aes(y=ratio, x=mean))+
#       geom_point(alpha=0.4)+
#       geom_hline(yintercept = 0, color="black")+
#       labs(title= paste0(sample_id,"-", data$library_layout[k]),
#            subtitle= paste0(study,"-","Old ref_ratio:",round(median(ase_df$ref_ratio),2)))
#     
#     p1=p0+
#       geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=ratio_q, group=1), color="red")+
#       geom_line(data=q_line[which(q_line$q=="low"),], aes(x=max, y=ratio_q, group=1), color="red")+
#       geom_line(data=positive_invert, aes(x=max, y=invert, group=1), color="purple")+
#       geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=median_ratio, group=1), color="green")+
#       geom_segment(data = line_segment, mapping = aes(x=min, y=1, xend=max, yend=1), inherit.aes = FALSE, color="blue",
#                    arrow = arrow(length = unit(0.1,"cm")))
#     print(p1)
#     p2=p0+
#       geom_quantile(quantiles = c(0.05, 0.5, 0.95),
#                     color="orange")
#     
#     print(p2)
#     
#     
#     p0=ggplot(ase_df, aes(y=log10(alt), x=log10(ref)))+
#       geom_point(alpha=0.4)+
#       geom_hline(yintercept = 0, color="black")+
#       labs(title= paste0(sample_id,"-", data$library_layout[k]),
#            subtitle= paste0(study,"-","Old ref_ratio:",round(median(ase_df$ref_ratio),2)))+
#       geom_quantile(quantiles = c(0.05, 0.5, 0.95),
#                     color="orange", )
#     
#     print(p0)
#     
#     
#     rm(ase_df)
#   }}
# dev.off()

#Train based on Geuvadis:
seq_mean=seq(0,4.6,by=0.1)
for(k in 1:nrow(geu)){
  sample_id<-geu$external_id[k]
  study<-geu$study[k]
  print(k)
  if(sample_id %in% bad_samples){ print("bad sample")}else{
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
    if(exists("ase_all")){
      ase_df<-ase_all
      rm(ase_all)
    } else{
      print("ase found")}
    ase_df<- ase_df %>% mutate(ratio=log10(alt) - log10(ref),
                               mean=(log10(alt) + log10(ref))/2)
    
    
    
    ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)
    
    
    
    q_line<-ase_df %>% group_by(mean_20tile) %>%
      reframe(median_ratio=median(ratio),
                enframe(quantile(ratio, c(0.05,0.95)), "quantile", "ratio_q"),
                min=min(mean),max=max(mean)) %>% 
      mutate(q=case_when(quantile== '5%' ~ "low",
                         quantile== '95%' ~ "high" ))
    q_line$sample<-sample_id
    
    if(k==1){
      q_line_geu<-q_line
    }else{
      q_line_geu<-rbind(q_line_geu,q_line)
    }

}}}


test_line<-q_line_geu %>%
  group_by(mean_20tile,q) %>%
  reframe(geu_mean=mean(ratio_q),
            enframe(quantile(ratio_q, c(0.05,0.95)), "CI", "CI_val")) %>%
  mutate(CI=case_when(CI== '5%' ~ "low_ci",
                     CI== '95%' ~ "high-ci" ))
test_line<-test_line[!is.na(test_line$mean_20tile),]
test_line$max<-sort(c(seq_mean[-c(1,2,3,4,5)],seq_mean[-c(1,2,3,4,5)],seq_mean[-c(1,2,3,4,5)],seq_mean[-c(1,2,3,4,5)]))
#saveRDS(test_line, "data/geuvadis_quantile.rds")

pdf(file="~/plot/ASE/test/quantile_adjust_geu-Test3.pdf", width = 10, height = 6) 
for(k in 1:20){
  sample_id<-data$external_id[k]
  study<-data$study[k]
  print(k)
  
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
    if(exists("ase_all")){
      ase_df<-ase_all
      rm(ase_all)
    } else{
      print("ase found")}
    ase_df<- ase_df %>% mutate(ratio=log10(alt) - log10(ref),
                               mean=(log10(alt) + log10(ref))/2)
    
  
    
    ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)

    q_line<-ase_df %>% group_by(mean_20tile) %>%
      reframe(median_ratio=median(ratio),
                enframe(quantile(ratio, c(0.05,0.95)), "quantile", "ratio_q"),
                min=min(mean),max=max(mean)) %>%
      mutate(q=case_when(quantile== '5%' ~ "low",
                         quantile== '95%' ~ "high" ))


    positive_invert<- q_line[which(q_line$q=="high"),] %>%
      mutate(dist= ratio_q - median_ratio,
             invert= median_ratio-dist)

    line_segment<-ase_df %>% group_by(mean_20tile) %>%
      summarise(min=min(mean),max=max(mean))

    test_line_sample<-left_join(test_line,q_line[,-6])
    test_line_sample$geu_mean[which(test_line_sample$q=="high")]<- test_line_sample$geu_mean[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
    test_line_sample$geu_mean[which(test_line_sample$q=="low")]<- test_line_sample$geu_mean[which(test_line_sample$q=="low")] + test_line_sample$median_ratio[which(test_line_sample$q=="low")]

    test_line_sample$CI_val[which(test_line_sample$q=="high")]<- test_line_sample$CI_val[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
    test_line_sample$CI_val[which(test_line_sample$q=="low")]<- test_line_sample$CI_val[which(test_line_sample$q=="low")] + test_line_sample$median_ratio[which(test_line_sample$q=="low")]

    p0=ggplot(ase_df, aes(y=ratio, x=mean))+
      geom_point(alpha=0.4)+
      geom_hline(yintercept = 0, color="black")+
      labs(title= paste0(sample_id,"-", data$library_layout[k]),
           subtitle= paste0(study,"-","Old ref_ratio:",round(median(ase_df$ref_ratio),2), "|red line is mean 95%/5% of geu samples, orange is 95%/5% mean"))

    p0=p0+
      geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=geu_mean, group=1), color="red")+
      geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=geu_mean, group=1), color="red")+
      geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
      geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
      geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
      geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
      geom_segment(data = line_segment, mapping = aes(x=min, y=1, xend=max, yend=1), inherit.aes = FALSE, color="blue",
                   arrow = arrow(length = unit(0.1,"cm")))
     #print(p1)
     # p2=p0+
     #   geom_quantile(quantiles = 0.95,
     #                 color="orange")+
     #   geom_quantile(quantiles = 0.05,
     #                 color="orange")
     #   # geom_quantile(quantiles = 0.5,
     #   #               color="green",method = "rqss", lambda = 0.1)
     # 
     # print(p2)


     p1=p0+
      geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=ratio_q, group=1), color="yellow")+
      geom_line(data=q_line[which(q_line$q=="low"),], aes(x=max, y=ratio_q, group=1), color="yellow")+
      #geom_line(data=positive_invert, aes(x=max, y=invert, group=1), color="salmon")+
      geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=median_ratio, group=1), color="green")
     #
    print(p1)

    
    rm(ase_df)
  }}
dev.off()

#---------------------------------------------------
#Test on Gtex
#---------------------------------------------------
geno_met<-read.csv("data/GTEx_geno_metadata.csv")
uni_norm<-NA
for(i in 10:length(geno_met$study)){
  study<-geno_met$study[3]
  print(study)
  
  ase_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rds"))
  
  ase_df<- ase_df %>% mutate(ratio=log10(alt) - log10(ref),
                             mean=(log10(alt) + log10(ref))/2)
  
  
  
  ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)
  
  q_line_all<-ase_df %>% group_by(mean_20tile,sample_id_rep) %>%
    reframe(median_ratio=median(ratio),
            enframe(quantile(ratio, c(0.05,0.95)), "quantile", "ratio_q"),
            min=min(mean),max=max(mean)) %>%
    mutate(q=case_when(quantile== '5%' ~ "low",
                       quantile== '95%' ~ "high" ))
  
  pdf(file=paste0("~/plot/ASE/test/",study,"222222.pdf"), width = 10, height = 6) 
  
  for(ss in 1:length(unique(q_line_all$sample_id_rep))){
    print(ss)
 q_line<-q_line_all[which(q_line_all$sample_id_rep %in% unique(q_line_all$sample_id_rep)[ss]),]
  test_line_sample<-left_join(test_line,q_line[,-7])
  test_line_sample<-test_line_sample %>% filter(q=="high",CI=="high-ci",
                                                max >=1,max <=2.5)
  
  test_line_sample$geu_mean[which(test_line_sample$q=="high")]<- test_line_sample$geu_mean[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
  test_line_sample$CI_val[which(test_line_sample$q=="high")]<- test_line_sample$CI_val[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
  
  
  test_line_sample<-test_line_sample %>% mutate(diff= ratio_q-CI_val) %>%  filter(diff>0)
  
  test_line_sample<-left_join(test_line,q_line[,-7])
  test_line_sample$geu_mean[which(test_line_sample$q=="high")]<- test_line_sample$geu_mean[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
  test_line_sample$geu_mean[which(test_line_sample$q=="low")]<- test_line_sample$geu_mean[which(test_line_sample$q=="low")] + test_line_sample$median_ratio[which(test_line_sample$q=="low")]
  
  test_line_sample$CI_val[which(test_line_sample$q=="high")]<- test_line_sample$CI_val[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
  test_line_sample$CI_val[which(test_line_sample$q=="low")]<- test_line_sample$CI_val[which(test_line_sample$q=="low")] + test_line_sample$median_ratio[which(test_line_sample$q=="low")]
  
  
  
  uni_norm[ss]<- mean(test_line_sample$diff,na.rm=T)
  
  plot<-ase_df %>% filter(sample_id_rep == unique(q_line_all$sample_id_rep)[ss])
  
  
  p0=ggplot(plot, aes(y=ratio, x=mean))+
    geom_point(alpha=0.4)+
    geom_hline(yintercept = 0, color="black")+
    labs(title= paste0(sample_id,"-", data$library_layout[k]),
         subtitle= paste0(study,"-","Old ref_ratio:",round(median(plot$ref_ratio),2), "Gtex, Thyroid"))
  
  p0=p0+geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=geu_mean, group=1), color="red")+
    geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=geu_mean, group=1), color="red")+
    geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
    geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
    geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
    geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")
    
  
  
  p1=p0+
    geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=ratio_q, group=1), color="yellow")+
    geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=median_ratio, group=1), color="green")+
    annotate("text", x = 3, y = 0.85, label = round(uni_norm[ss],3) )
  #
  print(p1)
  

}
dev.off()
which(uni_norm>0.12)




#---------------------------------------------------
data$uni_norm<-NA
for(k in 1:nrow(data)){
  sample_id<-data$external_id[k]
  study<-data$study[k]
  print(k)
  
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
    if(exists("ase_all")){
      ase_df<-ase_all
      rm(ase_all)
    } else{
      print("ase found")}
    ase_df<- ase_df %>% mutate(ratio=log10(alt) - log10(ref),
                               mean=(log10(alt) + log10(ref))/2)
    
    
    
    ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)
    
    q_line<-ase_df %>% group_by(mean_20tile) %>%
      reframe(median_ratio=median(ratio),
              enframe(quantile(ratio, c(0.05,0.95)), "quantile", "ratio_q"),
              min=min(mean),max=max(mean)) %>%
      mutate(q=case_when(quantile== '5%' ~ "low",
                         quantile== '95%' ~ "high" ))
    
    
    # positive_invert<- q_line[which(q_line$q=="high"),] %>%
    #   mutate(dist= ratio_q - median_ratio,
    #          invert= median_ratio-dist)
    # 
    # line_segment<-ase_df %>% group_by(mean_20tile) %>%
    #   summarise(min=min(mean),max=max(mean))
    
    test_line_sample<-left_join(test_line,q_line[,-6])
    test_line_sample<-test_line_sample %>% filter(q=="high",
                                                  CI=="high-ci",
                                                  max >=1,
                                                  max <=2.5)
    
    test_line_sample$geu_mean[which(test_line_sample$q=="high")]<- test_line_sample$geu_mean[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
    test_line_sample$CI_val[which(test_line_sample$q=="high")]<- test_line_sample$CI_val[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]

    
    test_line_sample<-test_line_sample %>% mutate(diff= ratio_q-CI_val) %>%  filter(diff>0)
    
    data$uni_norm[k]<- mean(test_line_sample$diff,na.rm=T)
    
    rm(ase_df)
  }}


data$uni_norm[1:20]
bad<-data[which(data$uni_norm>0.12),]
good<-data[-which(data$uni_norm>0.12),]
#----------------------------------------
#plot bad samples

pdf(file="~/plot/ASE/test/select_good.pdf", width = 10, height = 6) 
for(k in 1:20){
  sample_id<-good$external_id[k]
  study<-good$study[k]
  uni_norm<-round(good$uni_norm[k],3)
  print(k)
  
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))){
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
    if(exists("ase_all")){
      ase_df<-ase_all
      rm(ase_all)
    } else{
      print("ase found")}
    ase_df<- ase_df %>% mutate(ratio=log10(alt) - log10(ref),
                               mean=(log10(alt) + log10(ref))/2)
    
    
    
    ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)
    
    q_line<-ase_df %>% group_by(mean_20tile) %>%
      reframe(median_ratio=median(ratio),
              enframe(quantile(ratio, c(0.05,0.95)), "quantile", "ratio_q"),
              min=min(mean),max=max(mean)) %>%
      mutate(q=case_when(quantile== '5%' ~ "low",
                         quantile== '95%' ~ "high" ))
    
    
    positive_invert<- q_line[which(q_line$q=="high"),] %>%
      mutate(dist= ratio_q - median_ratio,
             invert= median_ratio-dist)
    
    line_segment<-ase_df %>% group_by(mean_20tile) %>%
      summarise(min=min(mean),max=max(mean))
    
    test_line_sample<-left_join(test_line,q_line[,-6])
    test_line_sample$geu_mean[which(test_line_sample$q=="high")]<- test_line_sample$geu_mean[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
    test_line_sample$geu_mean[which(test_line_sample$q=="low")]<- test_line_sample$geu_mean[which(test_line_sample$q=="low")] + test_line_sample$median_ratio[which(test_line_sample$q=="low")]
    
    test_line_sample$CI_val[which(test_line_sample$q=="high")]<- test_line_sample$CI_val[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
    test_line_sample$CI_val[which(test_line_sample$q=="low")]<- test_line_sample$CI_val[which(test_line_sample$q=="low")] + test_line_sample$median_ratio[which(test_line_sample$q=="low")]
    
    p0=ggplot(ase_df, aes(y=ratio, x=mean))+
      geom_point(alpha=0.4)+
      geom_hline(yintercept = 0, color="black")+
      labs(title= paste0(sample_id,"-", data$library_layout[k]),
           subtitle= paste0(study,"-","Old ref_ratio:",round(median(ase_df$ref_ratio),2), "|red line is mean 95%/5% of geu samples, orange is 95%/5% mean"))
    
    p0=p0+
      geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=geu_mean, group=1), color="red")+
      geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=geu_mean, group=1), color="red")+
      geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
      geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
      geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
      geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
      geom_segment(data = line_segment, mapping = aes(x=min, y=1, xend=max, yend=1), inherit.aes = FALSE, color="blue",
                   arrow = arrow(length = unit(0.1,"cm")))
    
    
    p1=p0+
      geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=ratio_q, group=1), color="yellow")+
      geom_line(data=q_line[which(q_line$q=="low"),], aes(x=max, y=ratio_q, group=1), color="yellow")+
      #geom_line(data=positive_invert, aes(x=max, y=invert, group=1), color="salmon")+
      geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=median_ratio, group=1), color="green")+
      annotate("text", x = 4, y = 1.25, label = uni_norm )
    #
    print(p1)
    
    
    rm(ase_df)
  }}
dev.off()

bad$uni_norm[which(bad$study=="DRP003703")]

