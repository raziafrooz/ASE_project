library(cowplot)

theme_set(theme_cowplot())

make_MA<-function(ase_df,test_line){
   seq_mean=seq(0,4.6,by=0.1)
    
    ase_df<- ase_df %>% mutate(ratio=log10(alt_count) - log10(ref_count),
                               mean=(log10(alt_count) + log10(ref_count))/2)
    
    
    
    ase_df$mean_20tile<-cut(ase_df$mean, seq_mean,include.lowest =T)
    
    q_line<-ase_df %>% group_by(mean_20tile) %>%
      reframe(median_ratio=median(ratio),
              enframe(quantile(ratio, c(0.05,0.95)), "quantile", "ratio_q"),
              min=min(mean),max=max(mean)) %>%
      mutate(q=case_when(quantile== '5%' ~ "low",
                         quantile== '95%' ~ "high" ))
    
    
    #positive_invert<- q_line[which(q_line$q=="high"),] %>%
    #  mutate(dist= ratio_q - median_ratio,
    #         invert= median_ratio-dist)
    
    # line_segment<-ase_df %>% group_by(mean_20tile) %>%
    #   summarise(min=min(mean),max=max(mean))
    
    test_line_sample<-left_join(test_line,q_line[,-6])
    test_line_sample$geu_mean[which(test_line_sample$q=="high")]<- test_line_sample$geu_mean[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
    test_line_sample$geu_mean[which(test_line_sample$q=="low")]<- test_line_sample$geu_mean[which(test_line_sample$q=="low")] + test_line_sample$median_ratio[which(test_line_sample$q=="low")]
    
    test_line_sample$CI_val[which(test_line_sample$q=="high")]<- test_line_sample$CI_val[which(test_line_sample$q=="high")] + test_line_sample$median_ratio[which(test_line_sample$q=="high")]
    test_line_sample$CI_val[which(test_line_sample$q=="low")]<- test_line_sample$CI_val[which(test_line_sample$q=="low")] + test_line_sample$median_ratio[which(test_line_sample$q=="low")]
    
    MA<-list(ase_df=ase_df,test_line_sample=test_line_sample,q_line=q_line)
    return(MA)
  }


plot_MA<-function(ase_df, sample_id,study,uni_norm,test_line_sample,q_line){
  p0=ggplot(ase_df, aes(y=ratio, x=mean))+
    geom_point(alpha=0.4)+
    geom_hline(yintercept = 0, color="black")+
    labs(title= sample_id,
         x="A", y="M")
         #subtitle= paste0(study,"-","Old ref_ratio:",round(median(ase_df$ref_ratio),2)))
  
  p0=p0+
    geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=geu_mean, group=1), color="red")+
    geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=geu_mean, group=1), color="red")+
    geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
    geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
    geom_line(data=test_line_sample[which(test_line_sample$q=="high"),][which(test_line_sample$CI=="low_ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")+
    geom_line(data=test_line_sample[which(test_line_sample$q=="low"),][which(test_line_sample$CI=="high-ci"),], aes(x=max, y=CI_val, group=1), color="orange", linetype="dashed")
  #geom_segment(data = line_segment, mapping = aes(x=min, y=1, xend=max, yend=1), inherit.aes = FALSE, color="blue",
  #            arrow = arrow(length = unit(0.1,"cm")))
  
  
  p1=p0+
    geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=ratio_q, group=1), color="magenta2")+
    geom_line(data=q_line[which(q_line$q=="low"),], aes(x=max, y=ratio_q, group=1), color="magenta2")+
    #geom_line(data=positive_invert, aes(x=max, y=invert, group=1), color="salmon")+
    geom_line(data=q_line[which(q_line$q=="high"),], aes(x=max, y=median_ratio, group=1), color="lightgreen")+
    annotate("text", x = 4, y = 1.25, label = paste0("regression_value:",round(uni_norm,2)) )
  #
  return(p1)
  
  
}
