#fwrite(gtex_w_ov, file="~/test/gtex_fixed.csv")
gtex_w_ov<-fread("~/test/gtex_fixed.csv")
wasp<-fread("/dcs07/hansen/data/recount_ASE/data/gtexVSrecount_wasp.csv.gz")
wasp_plot<- gtex_w_ov %>% mutate(sample=SAMPLE_ID) %>% select(sample,adj_ref_ratio) %>% right_join(wasp,by="sample")
wasp_plot$adj_ref_ratio[is.na(wasp_plot$adj_ref_ratio)]<-wasp_plot$recount3[is.na(wasp_plot$adj_ref_ratio)]
wasp_plot<-wasp_plot %>% pivot_longer(!c("sample","tissue"), names_to="pipeline",values_to="ref_ratio")


pdf(file="~/plot/ASE/adjustVSwasp.pdf", width = 10, height = 6)
ggplot(data=wasp_plot,aes(x=pipeline, y=ref_ratio))+
  geom_boxplot()+
  labs(title="After correcting overlap issue:14313 GTEx samples")
dev.off()

