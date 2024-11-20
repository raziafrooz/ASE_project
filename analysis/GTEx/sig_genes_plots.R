
library(tidyverse)
library(data.table)

recount_sig_snps<-fread("~/test/recount_sig_snps.csv.gz")
qc_df<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/gtex_qc_metadata.csv.gz")
qc_df<-qc_df %>% mutate(overlap= (star.average_input_read_length)-bc_frag.mode_length)


#colnames(recount_sig_snps)[1]<-"SAMPLE_ID"

qc_df$SAMPLE_ID<-str_sub(qc_df$external_id, end= -3)
qc_df<-qc_df %>% filter(SMGEBTCHT=="TruSeq.v1", overlap<200) %>%  dplyr::select(SAMPLE_ID,overlap )
try1<-right_join(qc_df,recount_sig_snps,by=c("SAMPLE_ID"))



pdf(file="~/plot/ASE/test_qc1.pdf", width = 10, height = 6)
ggplot(data=try1,aes(x=sig_gtex, y=overlap,color=SMNABTCHT))+
  geom_point(alpha=0.5)
dev.off()

pdf(file="~/plot/ASE/qc_plots1_2.pdf", width = 10, height = 6)
for(name_col in colnames(try1)[18:194]){
  
  print(name_col)
  
  pp=ggplot(try1,aes(y=.data[[name_col]], x=sig_gtex))+
    geom_point()+ 
    labs(title=name_col)
  
  print(pp)
}
dev.off()



pdf(file="~/plot/ASE/qc_plots2.pdf", width = 10, height = 6)

pp=ggplot(try1%>% filter(sig_gtex>3000),aes(x=SMGEBTCHD, y=sig_gtex))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust=1))+
  labs(title="There is a relationship with batch date and # of sig snps in gtex",
       subtitle="Here I am only showing samples with sig snps>3k")

print(pp)


class <- try1 %>%
  group_by(SMNABTCHT) %>% 
  summarize(label = paste0("n = ", n())) %>% ungroup()
pp=ggplot(try1,aes(x=SMNABTCHT, y=sig_gtex))+
  geom_jitter()+
  theme(axis.text.x = element_text(angle = 5, vjust = 0.5, hjust=1))+
  labs(title="There is a relationship with batch date and # of sig snps in gtex",
       subtitle=paste0("Total sample n=", nrow(try1)))+
  geom_text(data = class, aes(y = 30000, label = label))

print(pp)


pp=ggplot(try1,aes(color=SMNABTCHT, x=overlap,y=sig_gtex))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))+
  labs(title="There is a relationship with batch date and # of sig snps in gtex",
      subtitle=paste0("Total sample n=", nrow(try1)))

print(pp)

pp=ggplot(try1,aes(x=SMGEBTCHD, y=overlap))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 10, vjust = 0.5, hjust=1))+
  labs(title="There is a relationship with batch date and # of sig snps in gtex",
       subtitle=paste0("Total sample n=", nrow(try1)))

print(pp)
dev.off()



