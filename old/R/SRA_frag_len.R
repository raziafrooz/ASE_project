setwd("~/ASE/")
library(recount3)
library(tidyverse)
library(GenomicRanges)
library(ggplot2)

sra<-readRDS("~/plot/ASE/sra_qc.rds")
sra_norm<-readRDS("~/plot/ASE/sra_qc_normal.rds")

#id<-which(sra$ref_ratio<=0.5)
id<-which(sra$library_layout=="single")
s_id<-sample(id,500)
sra_subset<-sra[s_id,]

id<-which(sra$library_layout=="paired")
s_id<-sample(id,500)
sra_subset<-rbind(sra_subset,sra[s_id,])


s_id<-sample(nrow(sra_norm),500)
#sra_subset_norm<-sra_norm[s_id,]
sra_subset<-rbind(sra_subset,sra_norm[s_id,])


sra_subset$frag_len<-NA
for(ss in unique(sra_subset$study)){
  print(ss)
  url<-locate_url(
    ss,
    "data_sources/sra",
    type = "metadata")
  
  xx <-utils::read.delim(file_retrieve(url[4], verbose = FALSE))
  
  id<-match(sra_subset$sample_id,xx$external_id)
  id<-id[!is.na(id)]
  id2<-match(xx$external_id,sra_subset$sample_id)
  id2<-id2[!is.na(id2)]
  
  
  sra_subset$frag_len[id2]<-xx$bc_frag.mean_length[id]
  
if(sum(is.na(sra_subset$avg_len[id2]))!=0){
  
  xx <-utils::read.delim(file_retrieve(url[4], verbose = FALSE))

  id<-match(sra_subset$sample_id,xx$external_id)
  id<-id[!is.na(id)]
  id2<-match(xx$external_id,sra_subset$sample_id)
  id2<-id2[!is.na(id2)]
  
  
  sra_subset$avg_len[id2]<-xx$avg_len[id]
}
}



table(sra2$library_layout)
sra_subset$seq_overlap<- sra_subset$frag_len
indx<-sra_subset$library_layout=="paired"
sra_subset$seq_overlap[indx]<- sra_subset$frag_len[indx]-(sra_subset$avg_len[indx]*2)



plot_sra_subset<-sra_subset 
plot_sra_subset$fg_cut<-0
plot_sra_subset$fg_cut[plot_sra_subset$library_layout=="paired"]<- cut_number(plot_sra_subset$frag_len[plot_sra_subset$library_layout=="paired"],5)

library(ggridges)

pdf(file="~/plot/ASE/sra_frag_len.pdf", width = 10, height = 6)

ggplot(plot_sra_subset)+
  geom_violin(aes(x=factor(fg_cut),y=ref_ratio,color=library_layout))+
  geom_jitter(aes(x=factor(fg_cut),y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="random 10K SRA samples")

ggplot(plot_sra_subset)+
  geom_density_ridges(aes(frag_len,y=library_layout,fill=library_layout))+
  labs(title="random 10K SRA samples")

dev.off()


pdf(file="~/plot/ASE/SRA_library_layout.pdf", width = 10, height = 6)

ggplot(sra %>% filter(library_layout%in% c("paired","single")))+
  geom_violin(aes(x=library_layout,y=ref_ratio,color=library_layout))+
  geom_jitter(aes(x=library_layout,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)
dev.off()

#-----------------------------------------------

which(unique(sra_subset$study)=="SRP048664")
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
  
  
  
  lung<-sra_subset[id2,]
  xx$ref_ratio<-NA
  xx$ref_ratio[id]<-lung$ref_ratio
  xx<-xx[!is.na(xx$ref_ratio),]
  
  
  if(ss == unique(sra_subset$study)[1]){
    sra_qc<-xx}else{
      sra_qc<-rbind(sra_qc,xx)
    }
}

sra_qc$library_layout<-sra_subset$library_layout[match(sra_qc$external_id,sra_subset$sample_id)]
name_col <- colnames(sra_qc)[15]
pdf(file="~/plot/ASE/SRA_QC.pdf", width = 10, height = 6)
for(name_col in colnames(sra_qc)[4:55]){
  if( sum(duplicated(quantile(sra_qc[,name_col]))) == 0   ){
    
    print(name_col)
    plot2<- sra_qc %>% mutate(cut=cut_number(.data[[name_col]],6))
    
    
    
    p= ggplot(plot2)+
      geom_violin(aes(x=cut,y=ref_ratio, color=library_layout))+
      geom_jitter(aes(x=cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
      labs(title=name_col,
           x=paste0(name_col, "_cut"))
    print(p)
  }else{
    plot2<- sra_qc %>% mutate(ref_ratio_cut=cut_number(ref_ratio,6))
    
    p= ggplot(plot2)+
      geom_density_ridges(aes(x=.data[[name_col]],y=ref_ratio_cut, fill=library_layout), alpha=0.5)+
      labs(title=name_col)
    print(p)
 #   print("cannot make the plot")
  }
}
dev.off()

pdf(file="~/plot/ASE/SRA_QC3.pdf", width = 10, height = 6)
for(name_col in colnames(sra_qc)[109:113]){
  if( sum(duplicated(quantile(sra_qc[,name_col]))) == 0   ){
    
    print(name_col)
    plot2<- sra_qc %>% mutate(cut=cut_number(.data[[name_col]],6))
    
    
    
    p= ggplot(plot2)+
      geom_violin(aes(x=cut,y=ref_ratio, color=library_layout))+
      geom_jitter(aes(x=cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
      labs(title=name_col,
           x=paste0(name_col, "_cut"))
    print(p)
  }else{
    plot2<- sra_qc %>% mutate(ref_ratio_cut=cut_number(ref_ratio,6))
    
    p= ggplot(plot2)+
      geom_density_ridges(aes(x=.data[[name_col]],y=ref_ratio_cut, fill=library_layout), alpha=0.5)+
      labs(title=name_col)
    print(p)
    #   print("cannot make the plot")
  }
}
dev.off()
colnames(sra_qc)[17]
sra_qc$mean_mode<-sra_qc[,17]-sra_qc[,16]
sra_qc[900:920,16]
sra_qc[900:920,17]
sra_qc[900:920,16]-sra_qc[900:920,17]
sra_qc[,16]-sra_qc[,17]

sra_subset$frag_mode<-sra_qc$bc_frag.mode_length[match(sra_subset$sample_id,sra_qc$external_id)]
sra_subset$mode_overlap<- sra_subset$frag_mode
indx<-sra_subset$library_layout=="paired"
sra_subset$mode_overlap[indx]<- sra_subset$frag_mode[indx]-(sra_subset$avg_len[indx]*2)


pdf(file="~/plot/ASE/Kallisto_mode_frag.pdf", width = 10, height = 6)


plot2<-sra_qc %>% filter(library_layout=="paired") %>%  mutate( mean_mode_cut=cut_number(mean_mode,10),
                                                        ref_ratio_cut=cut_number(ref_ratio,10))
plot3<-sra_subset %>% filter(library_layout=="paired") %>%  mutate( mode_overlap_cut=cut_number(mode_overlap,10),
                                                                    seq_overlap_cut= cut_number(seq_overlap,10) )

p= ggplot(plot2)+
  geom_density_ridges(aes(x=mean_mode,y=ref_ratio_cut, fill=library_layout), alpha=0.5)+
  labs(title="using Kallisto:x axis is (mode.frag - mean.frag)")
print(p)


p= ggplot(plot3)+
  geom_violin(aes(x=mode_overlap_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=mode_overlap_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="Using mode frag to calculate (mode.frag - avg_len*2)")
print(p)

p= ggplot(plot3)+
  geom_violin(aes(x=seq_overlap_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=seq_overlap_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="Using mean frag to calculate (mean.frag - avg_len*2)")
print(p)


dev.off()


#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------

sra_subset_norm$library_layout<-NA
for(ss in unique(sra_subset_norm$study)){
  print(ss)
  url<-locate_url(
    ss,
    "data_sources/sra",
    type = "metadata")
  
  xx <-utils::read.delim(file_retrieve(url[1], verbose = FALSE))

  id<-match(sra_subset_norm$sample_id,xx$external_id)
  id<-id[!is.na(id)]
  id2<-match(xx$external_id,sra_subset_norm$sample_id)
  id2<-id2[!is.na(id2)]
  
  
  sra_subset_norm$library_layout[id2]<-xx$library_layout[id]
  
}



sra_subset_all<-rbind(sra_subset,sra_subset_norm)




for(ss in unique(sra_subset_norm$study)){
  print(ss)
  url<-locate_url(
    ss,
    "data_sources/sra",
    type = "metadata")
  
  xx <-utils::read.delim(file_retrieve(url[3], verbose = FALSE))
  
  id<-match(sra_subset_norm$sample_id,xx$external_id)
  id<-id[!is.na(id)]
  id2<-match(xx$external_id,sra_subset_norm$sample_id)
  id2<-id2[!is.na(id2)]
  
  
  
  lung<-sra_subset_norm[id2,]
  xx$ref_ratio<-NA
  xx$ref_ratio[id]<-lung$ref_ratio
  xx<-xx[!is.na(xx$ref_ratio),]
  
  
  if(ss == unique(sra_subset_norm$study)[1]){
    sra_qc_norm<-xx
    }else{
      sra_qc_norm<-rbind(sra_qc_norm,xx)
    }
}
sra_qc_norm[1:3,]
sra_qc_norm$library_layout<-sra_subset_norm$library_layout[match(sra_qc_norm$external_id,sra_subset_norm$sample_id)]




sra_subset_norm$frag_mode<-sra_qc_norm$bc_frag.mode_length[match(sra_subset_norm$sample_id,sra_qc_norm$external_id)]
sra_subset_norm$mode_overlap<- sra_subset_norm$frag_mode
indx<-which(sra_subset_norm$library_layout=="paired")
sra_subset_norm$mode_overlap[indx]<- sra_subset_norm$frag_mode[indx]-(sra_subset_norm$avg_len[indx]*2)


sra_subset_norm$seq_overlap<- sra_subset_norm$frag_len
indx<-which(sra_subset_norm$library_layout=="paired")
sra_subset_norm$seq_overlap[indx]<- sra_subset_norm$frag_len[indx]-(sra_subset_norm$avg_len[indx]*2)


pdf(file="~/plot/ASE/norm_Kallisto_mode_frag.pdf", width = 10, height = 6)


p= ggplot(sra_subset_norm)+
  geom_violin(aes(x=library_layout,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=library_layout,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="ref_ratio(0.48-0.52) vs library layout")
print(p)


plot3<-sra_subset_norm %>%  mutate( ref_ratio_cut=cut_number(ref_ratio,3))
p= ggplot(plot3)+
  geom_density_ridges(aes(x=frag_mode,y=ref_ratio_cut, fill=library_layout), alpha=0.5)+
  labs(title="ref_ratio(0.48-0.52) using Kallisto:x axis is mode.frag")
print(p)


plot3<-sra_subset_norm %>% filter(library_layout=="paired") %>%  mutate( mode_overlap_cut=cut_number(mode_overlap,10),
                                                                    seq_overlap_cut= cut_number(seq_overlap,10),
                                                                    ref_ratio_cut=cut_number(ref_ratio,4))





p= ggplot(plot3)+
  geom_violin(aes(x=mode_overlap_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=mode_overlap_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="Using mode frag to calculate (mode.frag - avg_len*2)")
print(p)

p= ggplot(plot3)+
  geom_violin(aes(x=seq_overlap_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=seq_overlap_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="Using mean frag to calculate (mean.frag - avg_len*2)")
print(p)


dev.off()


