sra_subset<-readRDS("~/plot/ASE/sra_subset.rds")
sra_subset<-sra_subset[,c(colnames(sra_subset)[1:4],"library_layout")]
colnames(sra_subset)

sra_subset$frag_len<-NA
sra_subset$seq_len<-NA
for(ss in unique(sra_subset$study)){
  print(ss)
  url<-locate_url(
    ss,
    "data_sources/sra",
    type = "metadata")
  
  # xx <-utils::read.delim(file_retrieve(url[3], verbose = FALSE))
  # 
  # id<-match(sra_subset$sample_id,xx$external_id)
  # id<-id[!is.na(id)]
  # id2<-match(xx$external_id,sra_subset$sample_id)
  # id2<-id2[!is.na(id2)]
  # 
  # 
  # sra_subset$frag_len[id2]<-xx$bc_frag.mode_length[id]
  # 
  # if(sum(is.na(sra_subset$seq_len[id2]))!=0){
  #   
    xx <-utils::read.delim(file_retrieve(url[4], verbose = FALSE))
    
    id<-match(sra_subset$sample_id,xx$external_id)
    id<-id[!is.na(id)]
    id2<-match(xx$external_id,sra_subset$sample_id)
    id2<-id2[!is.na(id2)]
    
    
    # sra_subset$seq_len[id2]<-xx$avg_len[id]
    sra_subset$num_bases[id2]<-xx$X.bases[id]
  }
}
#saveRDS(sra_subset,file="~/plot/test.rds")
sra_subset<-readRDS("~/plot/test.rds")
sra_subset$num_bases_mil<-sra_subset$num_bases/10^8
sra_subset[1:4,]
sra_subset$dis[which(is.na(sra_subset$dis))]<-"mid"
plot1<- sra_subset[!is.na(sra_subset$library_layout),] %>% filter(library_layout=="single")%>% mutate(seq_cut=cut_number(seq_len,4))
plot2<- sra_subset[!is.na(sra_subset$library_layout),] %>% filter(library_layout=="paired") %>% mutate(seq_cut=cut_number(seq_len,4),
                                                       frag_cut=cut_number(frag_len,4))

sra_subset$overlap<-sra_subset$frag_len-(2*sra_subset$seq_len)
plot3<- sra_subset[!is.na(sra_subset$library_layout),] %>% 
  group_by(library_layout) %>% mutate(seq_cut=cut_number(seq_len,4),overlap_cut=cut_number(overlap,4)) %>% 
  ungroup() %>% mutate(ref_ratio_cut=cut_number(ref_ratio,8))
table(plot3$ref_ratio_cut)
  
plot3<-plot3 %>% mutate(overlap_cut=cut_number(overlap,4))
table(plot3$overlap_cut,plot3$library_layout)
tail(sra_subset)
library(ggridges)
dev.off()

pdf(file="~/plot/ASE/SRA_overlap.pdf", width = 10, height = 6)
  
p= ggplot(plot3 %>% filter(library_layout=="paired")%>% mutate(overlap_cut=cut_number(overlap,4)))+
  geom_violin(aes(x=overlap_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=overlap_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="only looking at paired SRA samples: overlap is [frag_len-(2*seq_len)] ")
print(p)

p= ggplot(plot3)+
  geom_violin(aes(x=seq_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=seq_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="looking at 15k SRA samples: average seq length is shown")
print(p)

p= ggplot(plot3 %>% filter(library_layout=="paired")%>% mutate(frag_cut=cut_number(frag_len,4)))+
  geom_violin(aes(x=frag_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=frag_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="only looking at paired SRA samples: mod fragment length is shown ")
print(p)
dev.off()

p= ggplot(plot3)+
  geom_violin(aes(x=ref_ratio_cut,y=num_bases_mil, color=library_layout))+
  geom_jitter(aes(x=ref_ratio_cut,y=num_bases_mil,color=library_layout), alpha=0.5,width = 0.15)
print(p)

p= ggplot(plot3)+
  geom_density_ridges(aes(y=ref_ratio_cut,x=num_bases_mil, fill=library_layout), alpha=0.5)
print(p)

p= ggplot(plot3 %>% filter(library_layout=="paired")%>% mutate(overlap_cut=cut_number(overlap,4)))+
  geom_violin(aes(x=overlap_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=overlap_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)
print(p)

p= ggplot(plot3)+
  geom_density_ridges(aes(y=ref_ratio_cut,x=overlap, fill=library_layout), alpha=0.5)
print(p)

p= ggplot(plot3)+
  geom_density_ridges(aes(y=ref_ratio_cut,x=seq_len, fill=library_layout), alpha=0.5)
print(p)

p= ggplot(plot3)+
  geom_density_ridges(aes(y=ref_ratio_cut,x=frag_len, fill=library_layout), alpha=0.5)
print(p)

p= ggplot(plot3)+
  geom_violin(aes(x=dis,y=frag_len, color=library_layout))+
  geom_jitter(aes(x=dis,y=frag_len,color=library_layout), alpha=0.5,width = 0.15)
print(p)
  
p= ggplot(plot1)+
  geom_violin(aes(x=seq_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=seq_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)
print(p)

p= ggplot(plot1)+
  geom_point(aes(x=ref_ratio,y=seq_len, color=library_layout))
print(p)
    
    p= ggplot(plot2)+
      geom_violin(aes(x=seq_cut,y=ref_ratio, color=dis))+
      geom_jitter(aes(x=seq_cut,y=ref_ratio,color=dis), alpha=0.5,width = 0.15)

    print(p)
    p= ggplot(plot2)+
      geom_violin(aes(x=frag_cut,y=ref_ratio, color=dis))+
      geom_jitter(aes(x=frag_cut,y=ref_ratio,color=dis), alpha=0.5,width = 0.15)
    print(p)
    
    p= ggplot(plot3)+
      geom_point(aes(x=frag_len,y=seq_len, color=ref_ratio_cut,shape=library_layout),alpha=0.5)
    print(p)
    
  dev.off()
    
pdf(file="~/plot/ASE/SRA_seqVSfrag.pdf", width = 10, height = 6)
p= ggplot(sra_subset%>% filter(library_layout=="paired") %>% mutate(ref_ratio_cut=cut_number(ref_ratio,4)))+
  geom_point(aes(x=frag_len,y=seq_len,color=ref_ratio_cut,shape=library_layout),alpha=0.5)+
  labs(title="looking at paired SRA samples: average seq length vs mod fragment len")
print(p)

dev.off()

plot_df<- sra_subset%>% filter(library_layout=="single") %>% mutate(ref_ratio_cut=cut_number(ref_ratio,4),
                                                                    ref_ratio_cut2=ntile(ref_ratio,4))
table(plot_df$ref_ratio_cut2)
pdf(file="~/plot/ASE/SRA_seqlen.pdf", width = 10, height = 6)
p= ggplot(sra_subset[!is.na(sra_subset$library_layout),]%>% group_by(library_layout) %>% mutate(ref_ratio_cut=cut_number(ref_ratio,5)))+
  geom_density_ridges(aes(x=seq_len,y=ref_ratio_cut, fill=library_layout),alpha=0.5)+
  labs(title="looking at SRA samples: average seq length vs ref_ratio")
print(p)

p= ggplot(sra_subset[!is.na(sra_subset$library_layout),]%>% group_by(library_layout) %>% mutate(ref_ratio_cut=cut_number(ref_ratio,5)))+
  geom_density_ridges(aes(x=frag_len,y=ref_ratio_cut, fill=library_layout),alpha=0.5)+
  labs(title="looking at SRA samples: fragment length vs ref_ratio")
print(p)
dev.off()

plot_df<- sra_subset[!is.na(sra_subset$library_layout),] %>% 
  group_by(library_layout) %>% mutate(seq_cut=cut_number(seq_len,4),overlap_cut=cut_number(overlap,4),
                                      ref_ratio_cut=cut_number(ref_ratio,5)) %>% ungroup()
plot_df<- plot_df %>% filter(library_layout=="single")
table(plot_df$ref_ratio_cut,plot_df$seq_cut,plot_df$library_layout)

pdf(file="~/plot/ASE/SRA_overlap_v2.pdf", width = 10, height = 6)
  
p= ggplot(plot3 %>% filter(library_layout=="paired") %>% mutate(num_bases_cut= cut_number(num_bases_mil,6)))+
  geom_density_ridges(aes(y=num_bases_cut,x=overlap, fill=library_layout),alpha=0.5)+
  labs(title="only looking at paired SRA samples",
       subtitle="total number of bases across spots (not including both read mates if paired)(*10^8) vs overlap")
print(p)

p= ggplot(plot3 %>% filter(library_layout=="paired") %>% mutate(num_bases_cut= cut_number(num_bases_mil,6)))+
  geom_violin(aes(x=num_bases_cut,y=overlap, color=library_layout))+
  geom_jitter(aes(x=num_bases_cut,y=overlap,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="only looking at paired SRA samples",
       subtitle="total number of bases across spots (not including both read mates if paired)(*10^8) vs overlap")
print(p)


p= ggplot(plot3 %>% filter(library_layout=="paired") %>% mutate(num_bases_cut= cut_number(num_bases_mil,6)))+
  geom_violin(aes(x=num_bases_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=num_bases_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="only looking at paired SRA samples",
       subtitle="total number of bases across spots (not including both read mates if paired)(*10^8) vs overlap")
print(p)

ggplot(df)+
  geom_violin(aes(x=per_unmapped_tooShort_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=per_unmapped_tooShort_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="only looking at paired SRA samples",
       subtitle="%_of_reads_unmapped:_too_short vs overlap")

ggplot(df)+
  geom_violin(aes(x=per_unmapped_tooShort_cut,y=overlap, color=library_layout))+
  geom_jitter(aes(x=per_unmapped_tooShort_cut,y=overlap,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="only looking at paired SRA samples",
       subtitle="%_of_reads_unmapped:_too_short vs overlap")

dev.off()









p= ggplot(plot3 %>% mutate(num_bases_cut= cut_number(num_bases_mil,4)))+
  geom_violin(aes(x=num_bases_cut,y=ref_ratio, color=dis))+
  geom_jitter(aes(x=num_bases_cut,y=ref_ratio,color=dis), alpha=0.5,width = 0.15)
print(p)


p= ggplot(plot3 %>% mutate(num_bases_cut= cut_number(num_bases_mil,4)))+
  geom_violin(aes(x=num_bases_cut,y=overlap, color=dis))+
  geom_jitter(aes(x=num_bases_cut,y=overlap, color=dis), alpha=0.5,width = 0.15)
print(p)



p= ggplot(plot3 %>% mutate(num_bases_cut= cut_number(num_bases_mil,4)))+
  geom_violin(aes(x=num_bases_cut,y=overlap, color=library_layout))+
  geom_jitter(aes(x=num_bases_cut,y=overlap,color=library_layout), alpha=0.5,width = 0.15)
print(p)
dev.off()

#---------------------------------------------------------------------------------------
sra_subset$avg_input_read_len<- NA #star.average_input_read_len
sra_subset$unique_all<- NA #bc_auc.unique_reads_all_bases
sra_subset$per_unmapped_tooShort<-NA #star.._of_reads_unmapped._too_short

for(ss in unique(sra_subset$study)){
  print(ss)
  url<-locate_url(
    ss,
    "data_sources/sra",
    type = "metadata")
  
   xx <-utils::read.delim(file_retrieve(url[3], verbose = FALSE))
   
  id<-match(sra_subset$sample_id,xx$external_id)
  id<-id[!is.na(id)]
  id2<-match(xx$external_id,sra_subset$sample_id)
  id2<-id2[!is.na(id2)]


  sra_subset$per_unmapped_tooShort[id2]<-xx$star.._of_reads_unmapped._too_short[id]
  sra_subset$unique_all[id2]<-xx$bc_auc.unique_reads_all_bases[id]
  sra_subset$avg_input_read_len[id2]<-xx$star.average_input_read_length[id]
}

sra_subset$aFC<-abs(sra_subset$ref_ratio-0.5)



df<-sra_subset[!is.na(sra_subset$library_layout),] %>% group_by(library_layout) %>% mutate(aFC_cut=cut_number(aFC,10),
                          per_unmapped_tooShort_cut= cut_number(per_unmapped_tooShort,5),
                          unique_all_cut=cut_number(unique_all, 5),
                          avg_input_read_len_cut=cut_number(avg_input_read_len, 5),
                          ref_ratio_cut=cut_number(ref_ratio,5))

pdf(file="~/plot/ASE/SRA_test1.pdf", width = 10, height = 6)

ggplot(df)+
  geom_violin(aes(x=per_unmapped_tooShort_cut,y=aFC, color=library_layout))+
  geom_jitter(aes(x=per_unmapped_tooShort_cut,y=aFC,color=library_layout), alpha=0.5,width = 0.15)

ggplot(df)+
  geom_violin(aes(x=per_unmapped_tooShort_cut,y=overlap, color=library_layout))+
  geom_jitter(aes(x=per_unmapped_tooShort_cut,y=overlap,color=library_layout), alpha=0.5,width = 0.15)

ggplot(df)+
  geom_violin(aes(x=unique_all_cut,y=overlap, color=library_layout))+
  geom_jitter(aes(x=unique_all_cut,y=overlap,color=library_layout), alpha=0.5,width = 0.15)


ggplot(df)+
  geom_boxplot(aes(x=avg_input_read_len_cut,y=overlap, color=library_layout))+
  geom_jitter(aes(x=avg_input_read_len_cut,y=overlap,color=library_layout), alpha=0.5,width = 0.15)

dev.off()

pdf(file="~/plot/ASE/test.pdf", width = 10, height = 6)

ggplot(df)+
  geom_boxplot(aes(x=per_unmapped_tooShort_cut,y=seq_len, color=library_layout))+
  geom_jitter(aes(x=per_unmapped_tooShort_cut,y=seq_len,color=library_layout), alpha=0.5,width = 0.15)


p= ggplot(df)+
  geom_point(aes(x=frag_len,y=seq_len, color=ref_ratio_cut,shape=library_layout),alpha=0.5)
print(p)


p= ggplot(df)+
  geom_point(aes(x=unique_all,y=seq_len, color=ref_ratio_cut,shape=library_layout),alpha=0.5)
print(p)

p= ggplot(df)+
  geom_point(aes(x=per_unmapped_tooShort,y=seq_len, color=ref_ratio_cut,shape=library_layout),alpha=0.5)
print(p)

p= ggplot(df)+
  geom_point(aes(x=avg_input_read_len,y=seq_len, color=ref_ratio_cut,shape=library_layout),alpha=0.5)
print(p)
dev.off()


