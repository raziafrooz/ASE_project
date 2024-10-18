library(data.table)
library(tidyverse)
library(MetBrewer)

cc_fill <-scale_fill_manual(values=c(met.brewer("Redon"),met.brewer("Isfahan1"),met.brewer("Isfahan2")))
cc_color <-scale_color_manual(values=c(met.brewer("Redon"),met.brewer("Isfahan1"),met.brewer("Isfahan2")))

ase_df<-fread("/dcs07/hansen/data/recount_ASE/data/sra_ASE_stat.csv.gz")
recount3_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/metadata/Recount3_metadata.tsv", header= T, sep = "\t",quote="")



ase_df$library_layout<-recount3_metadata$library_layout[match(ase_df$experiment_acc,recount3_metadata$experiment_acc)]
ase_df$bc_frag.mode_length<-recount3_metadata$bc_frag.mode_length[match(ase_df$experiment_acc,recount3_metadata$experiment_acc)]
#ase_df$avg_len<-as.numeric(recount3_metadata$avg_len[match(ase_df$experiment_acc,recount3_metadata$experiment_acc)])
ase_df$star.average_input_read_length<-as.numeric(recount3_metadata$star.average_input_read_length[match(ase_df$experiment_acc,recount3_metadata$experiment_acc)])

ase_df<-ase_df[-which(ase_df$ref_ratio==0),]



ase_df<-ase_df %>% mutate(overlap= star.average_input_read_length-bc_frag.mode_length)

#Single end should not have bc_frag.mode_length, but some have >0 values--> fix them
ase_df$overlap[which(ase_df$library_layout=="single" & ase_df$bc_frag.mode_length==0)]<-ase_df$star.average_input_read_length[which(ase_df$library_layout=="single" & ase_df$bc_frag.mode_length==0)]

#-----------------------------
#Remove samples that are not bulk and added by mistake and also
#remove cancer samples for now


scRNA<-readRDS("/dcs07/hansen/data/recount_ASE/data/potential_single_cell.rds")
cancer_df<-readRDS("/dcs07/hansen/data/recount_ASE/data/cancer_annot.rds")
cancer_df$external_id<-recount3_metadata$external_id[match(cancer_df[,1],recount3_metadata$sample_acc)]

bad_samp<-data.frame(external_id=unique(c(scRNA,cancer_df$external_id)))
dim(bad_samp) #106,182
bad_samp$experiment_acc<-recount3_metadata$experiment_acc[match(bad_samp$external_id,recount3_metadata$external_id)]


ase_df<-ase_df[-which(ase_df$experiment_acc %in% unique(bad_samp$experiment_acc)),]
#-----------------------------




#--------------------------------------------
#SRA ref_ratio distribution plot:
#--------------------------------------------
pdf(file="~/plot/ASE/ref_ratio_sra.pdf", width = 6, height = 4.5)
ggplot(ase_df,aes(x=library_layout, y= ref_ratio, color=library_layout, fill=library_layout))+
  geom_jitter(size=1, alpha=0.5)+
  cc_fill+
  cc_color
  
ggplot(ase_df,aes(x=ref_ratio, color=library_layout, fill=library_layout))+
  geom_histogram(alpha=0.5)+
  cc_fill+
  cc_color

ggplot(ase_df[which(ase_df$library_layout!="single"),],aes(x=overlap, color=library_layout, fill=library_layout))+
  geom_histogram(alpha=0.5)+
  cc_fill+
  cc_color
dev.off()

#--------------------------------------------

single_df <-ase_df %>% filter(library_layout=="single") %>% 
  mutate(overlap_ntile=as.factor(cut_number(overlap,5)),
         avg_len_ntile= as.factor(cut_number(star.average_input_read_length,5)))

paired_df<-ase_df %>% filter(library_layout!="single") %>% 
  mutate(overlap_ntile=as.factor(cut_number(overlap,10)),
         aFC= abs(0.5 - ref_ratio))



pdf(file="~/plot/ASE/overlap_sra2.pdf", width = 8, height = 4.5)

p1= ggplot(paired_df)+
    cc_fill+
    cc_color+
  labs(title="Paired-end samples:\noverlap is star.average_input_read_length-bc_frag.mode_length",
       subtitle=paste0("There are ",nrow(paired_df), " pared_end seq samples"))

p2=p1+
  geom_violin(aes(x=overlap_ntile, y= ref_ratio, color=overlap_ntile, fill=overlap_ntile))+
  ylim(c(0.3,0.7))
  print(p2)

  
  p3= p1+
    geom_point(aes(x=star.average_input_read_length, y= ref_ratio, color=overlap_ntile, fill=overlap_ntile), alpha=0.4)+
    xlim(c(0,350))
  
  print(p3)


p3= p1+
  geom_point(aes(x=overlap, y= ref_ratio, color=overlap_ntile, fill=overlap_ntile), alpha=0.4)+
  xlim(c(-22,100))

print(p3)



p_single= ggplot(single_df)+
  cc_fill+
  cc_color+
  labs(title="Single-end samples: average length vs ref_ratio",
       subtitle=paste0("There are ",nrow(single_df), " single_end seq samples"))

p1.1=p_single+
  geom_violin(aes(x=avg_len_ntile, y= ref_ratio, color=avg_len_ntile, fill=avg_len_ntile))+
  ylim(c(0.45,0.6))
 
print(p1.1) 

p1.2=p_single+ 
  geom_point(aes(x=star.average_input_read_length, y= ref_ratio, color=avg_len_ntile, fill=avg_len_ntile),size=1, alpha=0.5)
print(p1.2) 

dev.off()



#------------------------------------------------------------------------------------------------------------------
recount3_metadata$ref_ratio<-as.numeric(ase_df$ref_ratio[match(recount3_metadata$experiment_acc,ase_df$experiment_acc)])


n<-sample(nrow(recount3_metadata), 20000)
plot_df<-recount3_metadata[n,]
colnames(plot_df)[40:50]
pdf(file="~/plot/ASE/sra_qc_plots.pdf", width = 10, height = 6)
for(name_col in colnames(plot_df)[47:198]){
  
  print(name_col)
  
  pp=ggplot(plot_df,aes(x=.data[[name_col]], y= ref_ratio, color = study))+
    geom_point()+ theme(legend.position="none")+
    labs(title=name_col)
  
  print(pp)
}
dev.off()

#------------------------------------------------------------------------------------------------------------------
metadata<-read.csv("/dcs07/hansen/data/recount_ASE/metadata/ASE_metadata.csv")



paired<-ase_df %>%
  filter(!is.na(ref_ratio),library_layout!="single") %>%
  mutate(ratio_ntile=ntile(ref_ratio,5)) %>% 
  group_by(ratio_ntile) %>% 
  sample_n(10)

pdf(file="~/plot/ASE/sra_indv_ase2.pdf", width = 10, height = 6)

for (i in 1:nrow(paired)){
  print(i)
  exp_id <- paired$experiment_acc[i]
  
  ase_indv<-fread(metadata$ASE_path[metadata$experiment_acc==exp_id])
  

  
  pp=ggplot(ase_indv,aes(x=log2(ref_count), y= log2(alt_count)))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1, color="red")+
    labs(title=paste0("SRA counts for ",exp_id,", nSNP= ", nrow(ase_indv)),
         subtitle=paste0("ref_ratio = " ,round(median(ase_indv$ref_ratio),3) ))
  
  print(pp)
  
}
  
  dev.off()
  
#-----------------------------------------------------------
  quantile(paired$overlap, probs = seq(0, 1, by = .1))
  quantile(single$overlap, probs = seq(0, 1, by = .1))
  paired<-ase_df %>%
    filter(!is.na(ref_ratio),library_layout!="single") %>%
    mutate(ratio_ntile=ntile(overlap,5)) %>% 
    group_by(ratio_ntile) %>% 
    sample_n(10)
  
  
  paired<-ase_df %>%
    filter(!is.na(ref_ratio),library_layout!="single", overlap>200)%>% 
    sample_n(30)
  
  
  pdf(file="~/plot/ASE/test4.pdf", width = 10, height = 6)
  
  for (i in 1:nrow(paired)){
    print(i)
    exp_id <- paired$experiment_acc[i]
    
    ase_indv<-fread(metadata$ASE_path[metadata$experiment_acc==exp_id])
    
    ov<-paired$overlap[i]
    
    pp=ggplot(ase_indv,aes(x=log2(ref_count), y= log2(alt_count)))+
      geom_point()+
      geom_abline(intercept = 0, slope = 1, color="red")+
      labs(title=paste0("SRA counts for ",exp_id,", nSNP= ", nrow(ase_indv)),
           subtitle=paste0("ref_ratio = " ,round(median(ase_indv$ref_ratio),3),
                           ", overlap = ", ov))
    
    print(pp)
    
  }
  
  dev.off()
  
  
  