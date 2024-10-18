geno_met<-read.csv("data/GTEx_geno_metadata.csv")
for(i in 1:length(geno_met$study)){
  study<-geno_met$study[i]
  print(study)
  study_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rds") )
  
  df<-study_df %>% filter(true_genotype!=2) %>% select(chr,start)
  if(i==1){
    blacklist<-df
  }else{
    blacklist<-rbind(blacklist,df)
  }
}  


blacklist_final <-blacklist %>% group_by(chr, start) %>% summarize(n=n()) %>% filter(n>20) %>% ungroup()
saveRDS(blacklist_final, file= "data/blacklist.rds")

indx<-paste0(ase_df$chr,ase_df$pos)
indx2<-paste0(blacklist_final$chr,blacklist_final$start)



pdf(file="~/plot/ASE/test.pdf", width = 10, height = 6)
ggplot(ase_df)+
  geom_point(aes(y=log2(alt),x=log2(ref)), alpha=0.5)

ggplot(ase_df[-which(indx%in%indx2),])+
  geom_point(aes(y=log2(alt),x=log2(ref)), alpha=0.2)+
  geom_point(data=ase_df[which(indx%in%indx2),],aes(y=log2(alt),x=log2(ref)), alpha=0.2, color="red")


ggplot(ase_df[which(indx%in%indx2),])+
  geom_point(aes(y=log2(alt),x=log2(ref)))

ggplot(ase_df[which(indx%in%indx2),])+
  geom_histogram(aes(ref_ratio))

ggplot(ase_df[-which(indx%in%indx2),])+
  geom_histogram(aes(ref_ratio))

ggplot(ase_df)+
  geom_histogram(aes(ref_ratio))

dev.off()


