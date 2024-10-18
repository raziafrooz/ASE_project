sra_subset<-readRDS("~/plot/ASE/sra_subset.rds")
sra_subset[1,]
sra_met<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/metadata/all_SRA.csv")
sra_met[1,]
recount3_chr_mapping <- "/dcl02/lieber/ajaffe/recount-pump/recount3.alts.chromosome_mappings.tsv"
snps_gr<-readRDS("/dcs04/hansen/data/recount_genotype/biallelic_snps/biallelic_SNP_gr.rds")
sra_geno<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/all_SRA.csv")

k=131241
which(sra_met$sample_id==sra_subset$sample_id[1])
geno2 %>% group_by(pred_genotype) %>% summarize(n())
geno3
coverage_cutoff=8
for(k in 1067:nrow(sra_subset)){
  print(k)
  mm<-which(sra_met$sample_id==sra_subset$sample_id[k])
  sample_id<-sra_met$sample_id[mm]
  study<-sra_met$study[mm]
  bigWig_path <- sra_met$total[mm]
  alt_path<-sra_met$alt[mm]
    # geno<-as_tibble(read.csv(sra_geno$genotypedSamples[k]))
    # geno <-geno %>% filter(pred_genotype==2, coverage>=8)
    # if(nrow(geno)>0){
    #   
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id,"_ase.rda"))
    if(exists("ase_df"))
    {
    geno <-ase_df#readRDS(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, ".rds"))
    rm(ase_df)
    }else{
      geno <-ase_all
      rm(ase_all)
      }
      geno_gr<-makeGRangesFromDataFrame(geno,seqnames="chr",start.field ="pos",end.field = "pos")
      genome(geno_gr)<-"hg38"
      ov<-findOverlaps(geno_gr,snps_gr )
      snp_gr<-snps_gr[subjectHits(ov)]      
      
      
      
      #----------------------------
      #read alt
      
      
      cat("Loading in: ", alt_path, "\n")
      temp_folder="~/test"
      temp_altFile <- paste0(temp_folder, sample_id, ".alt.temp.csv")
      tryCatch({
        system(paste0("zstdcat ", alt_path, " > ", temp_altFile))
      },
      error = function(e){
        message("Error reading alt file:\n", e)
        message("'zstdcat' software is not installed")
      },
      warning = function(w){
        message("Warning reading alt file:\n", w)
        message("'zstdcat' software is not installed")
      },
      finally = {
      })
      
      
      alt <- data.table::fread(temp_altFile, select = c(1, 2, 4))
      system(paste0("rm ", temp_altFile))
      if(nrow(alt) == 0) {
        return(NA)
      }
      colnames(alt) <- c("chr", "pos", "alt")
      #collapse alt counts.
      alt_count <- alt[, .N, by = .(chr, pos, alt)]
      rm(alt)
      alt_count <- alt_count[!is.na(alt_count$alt) ,]
      #`alt_count` is 0-based, so we add 1 to positions.
      alt_count$pos <- alt_count$pos + 1
      #use chromosome names that are used in recount3.
      chr_mapping <- fread(recount3_chr_mapping, header = FALSE)
      alt_count$chr <- chr_mapping$V2[match(alt_count$chr, chr_mapping$V1)]
      alt_count_gr <- GenomicRanges::GRanges(seqnames = alt_count$chr,
                                             ranges = IRanges(alt_count$pos,
                                                              alt_count$pos))
      ov <- GenomicRanges::findOverlaps(alt_count_gr, snp_gr)
      #Subset `alt_count_gr`, `alt_count` to the positions from SNPs, and compute
      #its own unique key.
      alt_count_gr <- alt_count_gr[queryHits(ov)]
      alt_count <- alt_count[queryHits(ov)]
      alt_key <- paste(as.character(seqnames(alt_count_gr)), start(alt_count_gr),
                        snp_gr$ref_seq[subjectHits(ov)], alt_count$alt, sep = "_")
      # #Even though `alt_key` is in the same positions as `filtered_snps_key`, many of the
      # #`alt_key` entries refer to alternate alleles that we are not tracking.
      # #We need to match a second time, this time using the entire key.
      filtered_snps_key <- paste(as.character(seqnames(snp_gr)), start(snp_gr),
                                 snp_gr$ref_seq, snp_gr$alt_seq, sep = "_")
      idx <- match(filtered_snps_key, alt_key)
      snps_key_idx <- which(!is.na(idx))
      alt_key_idx <- idx[!is.na(idx)]

      # final_alt_count <- rep(0, length(snp_gr))
      # final_alt_count[snps_key_idx] <- alt_count[alt_key_idx]$N
      # 
      #Compute M and S values.
      #ref_count <- bigwig_count - final_alt_count
      #We sometimes have
      #all.equal(geno$alt, final_alt_count)
      
      alt_count_gr_err<-alt_count_gr[-alt_key_idx,]
      alt_count_err<-alt_count[-alt_key_idx,]
      
      dd <- findOverlaps(alt_count_gr_err,snp_gr)
      dd_2<-dd[which(duplicated(subjectHits(dd)))]
      dd_3<-dd_2[which(duplicated(subjectHits(dd_2)))]
      
      error1<-dd[which(!duplicated(subjectHits(dd)))]
      error2<-dd_2[which(!duplicated(subjectHits(dd_2)))]
      error3<-dd_3[which(!duplicated(subjectHits(dd_3)))]
      
      geno$alt_err1<-NA
      geno$alt_err1_count<-0
      
      geno$alt_err2<-NA
      geno$alt_err2_count<-0
      
      geno$alt_err3<-NA
      geno$alt_err3_count<-0
      
      geno$alt_err1[subjectHits(error1)]<-alt_count_err$alt[queryHits(error1)]
      geno$alt_err1_count[subjectHits(error1)]<-alt_count_err$N[queryHits(error1)]
      
      geno$alt_err2[subjectHits(error2)]<-alt_count_err$alt[queryHits(error2)]
      geno$alt_err2_count[subjectHits(error2)]<-alt_count_err$N[queryHits(error2)]
      
      geno$alt_err3[subjectHits(error3)]<-alt_count_err$alt[queryHits(error3)]
      geno$alt_err3_count[subjectHits(error3)]<-alt_count_err$N[queryHits(error3)]
      
      

      saveRDS(geno,paste0("~/plot/test/", sample_id,"_ase_error.rds"))
}

sra_subset$nHet<-NA
sra_subset$nError<-NA
sra_subset$sum_error<- NA
sra_subset$sum_inser<-NA
sra_subset$total_cov<-NA
for(k in 1:1067){
  print(k)
sample_id<-sra_subset$sample_id[k]


geno<-readRDS(paste0("~/plot/test/", sample_id,"_ase_error.rds"))
sra_subset$nHet[k]<-nrow(geno)
sra_subset$nError[k]<-sum(geno$alt_err1_count>0)

sra_subset$sum_error[k] <- sum(geno$alt_err1_count+geno$alt_err2_count+geno$alt_err3_count)
sra_subset$total_cov[k] <- sum(geno$total)
sra_subset$sum_inser[k] <- length(which(nchar(geno$alt_err1)>1))

}
sra_subset$perc_error<-(sra_subset$sum_error/sra_subset$total_cov)*100
sra_subset$aFC<-abs(sra_subset$ref_ratio-0.5)


#-----------------------------------------------
#Recalc total and ref counts

recalc<-data.frame(old=sra_subset$ref_ratio, new=rep(NA, nrow(sra_subset)))
k=1
for(k in 1:1067){
  print(k)
  sample_id<-sra_subset$sample_id[k]
  
  
  geno<-readRDS(paste0("~/plot/test/", sample_id,"_ase_error.rds"))
  
  geno<-geno %>%
    rowwise() %>%
    mutate(alt= round(alt,0),
             new_total=total-sum(alt_err1_count,alt_err2_count,alt_err3_count)) %>% 
    mutate(new_ref=new_total-alt ) %>% 
    filter(new_total>=8, new_ref>=0) %>% 
    mutate(new_ref_ratio=new_ref/new_total)
  
  recalc$new[k]<-median(geno$new_ref_ratio)
  
  
}


pdf(file="~/plot/ASE/recalc_miscount.pdf", width = 10, height = 6)

p= ggplot(recalc)+
  geom_point(aes(x=old,y=new),alpha=0.5)+
  geom_abline(intercept = 0, slope = 1, color="red",linetype = "dashed",alpha=0.3)+
  labs(title="~1k SRA samples randomly selected (mostly extreme ref_ratio)",
       subtitle= "new ref ratio was calc after substracting total-error and recalc ref counts")
print(p)
dev.off()


ase_all_gr<-makeGRangesFromDataFrame(ase_all, seqnames.field = "chr", start.field = "pos", end.field = "pos")
ov<-findOverlaps(bed_graph,ase_all_gr)

x<-ase_all[unique(subjectHits(ov)),]

#-----------------------------------------------
#mappability
#-----------------------------------------------
recalc$after_mappability<-NA
#I found a bed file including the problematic sites from ENCODE
#https://github.com/Boyle-Lab/Blacklist/tree/master
black_list<-import("~/ASE/data/hg38-blacklist.v2.bed.gz")  

#https://bismap.hoffmanlab.org
bed_graph<-import("~/k50.Unique.Mappability.bb")

for(k in 1:1067){
  print(k)
  sample_id<-sra_subset$sample_id[k]
  
  
  geno<-readRDS(paste0("~/plot/test/", sample_id,"_ase_error.rds"))

  
  geno<-geno %>%
    rowwise() %>%
    mutate(alt= round(alt,0),
           new_total=total-sum(alt_err1_count,alt_err2_count,alt_err3_count)) %>% 
    mutate(new_ref=new_total-alt ) %>% 
    filter(new_total>=8, new_ref>=0) %>% 
    mutate(new_ref_ratio=new_ref/new_total)
  
 
  
  
  geno_gr<-makeGRangesFromDataFrame(geno,seqnames="chr",start.field ="pos",end.field = "pos")
  genome(geno_gr)<-"hg38"
  ov<-findOverlaps(bed_graph,geno_gr)
  
  geno<-geno[unique(subjectHits(ov)),]
  
  geno_gr<-makeGRangesFromDataFrame(geno,seqnames="chr",start.field ="pos",end.field = "pos")
  genome(geno_gr)<-"hg38"
  ov<-findOverlaps(black_list,geno_gr)
  
  geno<-geno[-unique(subjectHits(ov)),]
  
  recalc$after_mappability[k]<-median(geno$new_ref_ratio)
  
}
recalc$perc_error<-sra_subset$perc_error

recalc[1:10,]



#-----------------------------------------------
#plot few samples
#-----------------------------------------------


c(sample(1:1067,10 ))

k=776

pdf(file="~/plot/ASE/miscount_ind_sample.pdf", width = 10, height = 6)
for(k in c(58,33,638,108,344,1010,237,697,923,68)){
  print(k)
  sample_id<-sra_subset$sample_id[k]
  perc_error<-sra_subset$perc_error[k]
  
  geno<-readRDS(paste0("~/plot/test/", sample_id,"_ase_error.rds"))
  
  geno<-geno %>%
    rowwise() %>%
    mutate(alt= round(alt,0),
           new_total=total-sum(alt_err1_count,alt_err2_count,alt_err3_count)) %>% 
    mutate(new_ref=new_total-alt,
           new_ref_ratio= new_ref/ new_total) 
  
new_ref_ratio<- round( median(geno$new_ref_ratio[-which(geno$new_ref<0 | geno$new_total<8)]) ,3)
log2(geno$ref_ratio[1:4])+log2(geno$new_ref_ratio[1:4])/2
  x<-geno %>%  
    summarize(ratio=(log2(ref_ratio)-log2(new_ref_ratio)),
           mean=((log2(ref_ratio)+log2(new_ref_ratio))/2 )
           )
  
  text<-c(paste0("percent error=",round(perc_error,2)),
          paste0("new_ref_IsNeg=",sum(geno$new_ref<0) ),
          paste0("new_total_IsLess8=",sum(geno$new_total<8)),
          paste0("median_new=",new_ref_ratio ),
          paste0("median_old=",round(median(geno$ref_ratio),3))
  )
  
  p= ggplot(x,aes(x=mean,y=ratio))+
    geom_point(alpha=0.5)+
    geom_smooth()+
    labs(title= sample_id,
         subtitle="MA plot for old ref_ratio vs new (percent_error=(colsum(error)/total))")+
    annotate("text", x = -1, y = c(0.2,0.35,0.45,0.55,0.65), label = text)
  print(p)
  
  p= ggplot(geno,aes(x=ref_ratio,y=new_ref_ratio))+
    geom_point(alpha=0.5)+
    labs(title= sample_id,
         subtitle="plot for old ref_ratio vs new (percent_error=(colsum(error)/total))")+
    annotate("text", x = 0.4, y = c(0.5,0.525,0.55,0.575,0.6), label = text)+
    geom_abline(intercept = 0, slope = 1, color="red",linetype = "dashed",alpha=0.3)
  print(p)

  
}

pdf(file="~/plot/ASE/test1.pdf", width = 10, height = 6)
p= ggplot(geno, aes(new_ref_ratio))+
  geom_histogram( alpha=0.5)+
  xlim(c(0,1))
print(p)

dev.off()

sra_subset$ntile<-ntile(sra_subset$perc_error,10)
table(sra_subset$library_layout, sra_subset$ntile,sra_subset$interval)




total_dif<-data.frame(old_total=rep(NA, nrow(sra_subset)), 
                      new_total=rep(NA, nrow(sra_subset)),
                      perc_error=rep(NA, nrow(sra_subset)))
k=669
for(k in 669:1067){
  print(k)
  sample_id<-sra_subset$sample_id[k]
  
  
  geno<-readRDS(paste0("~/plot/test/", sample_id,"_ase_error.rds"))
  
  geno<-geno %>%
    rowwise() %>%
    mutate(alt= round(alt,0),
           error=sum(alt_err1_count,alt_err2_count,alt_err3_count),
           new_total=total-error) #%>% 
    # mutate(new_ref=new_total-alt ) %>% 
    # filter(new_total>=8, new_ref>=0) %>% 
    # mutate(new_ref_ratio=new_ref/new_total)
  
  total_dif$old_total[k]<-sum(geno$total)
  total_dif$new_total[k]<-sum(geno$new_total)
  total_dif$perc_error[k]<-sum(geno$error)
  print(sum(geno$new_total<0))
}

total_dif$total_diff<-(total_dif$perc_error/total_dif$old_total)*100
total_dif$ref_ratio<-sra_subset$ref_ratio
total_dif[1,]

pdf(file="~/plot/ASE/total_diff.pdf", width = 10, height = 6)
p= ggplot(total_dif)+
  geom_point(aes(x=total_diff,y=perc_error), alpha=0.5)+
  labs(x="percent error", y="difference between new total and old total",
       title="1K SRA sample")
print(p)

p= ggplot(total_dif)+
  geom_point(aes(x=total_diff,y=perc_error), alpha=0.5)+
  xlim(c(0,2))+
  ylim(c(0,10000))+
  labs(x="percent error", y="difference between new total and old total",
       title="1K SRA sample (no outliers)")
print(p)

p= ggplot(total_dif)+
  geom_point(aes(x=total_diff,y=log2(old_total)), alpha=0.5)+
  labs(x="percent error", y="log2(old total)",
       title="1K SRA sample")+
  xlim(0,5)
print(p)

p= ggplot(total_dif,aes(x=total_diff,y=abs(ref_ratio-0.5)))+
  geom_point(alpha=0.5)+
  labs(x="percent error", y="aFC",
       title="1K SRA sample")+
  geom_smooth()
print(p)

p= ggplot(total_dif,aes(x=total_diff,y=abs(ref_ratio-0.5)))+
  geom_point(alpha=0.5)+
  labs(x="percent error", y="aFC",
       title="1K SRA sample (no outlier)")+
  xlim(0,1)+
  geom_smooth()
print(p)

p= ggplot(total_dif,aes(x=total_diff,y=ref_ratio))+
  geom_point(alpha=0.5)+
  labs(x="percent error", y="ref_ratio",
       title="1K SRA sample (no outlier)")+
  xlim(0,1)+
  geom_smooth()
print(p)
dev.off()


#-----------------------------------------------
#plot
#-----------------------------------------------




plot<-sra_subset %>%  filter(library_layout=="paired") %>% 
  mutate(ref_ratio_cut=cut_number(ref_ratio,4),
         aFC_cut=ntile(aFC,4),
         nHet_cut=nError/nHet,
         perc_error_cut=cut_number(perc_error,10))
library(ggridges)


pdf(file="~/plot/ASE/test.pdf", width = 10, height = 6)
p= ggplot(plot)+
  geom_point(aes(x=nHet_cut,y=ref_ratio, color=library_layout), alpha=0.5)
print(p)

p= ggplot(plot)+
  geom_point(aes(x=nHet,y=ref_ratio, color=library_layout), alpha=0.5)
print(p)


p= ggplot(plot)+
  geom_point(aes(x=nError,y=ref_ratio, color=library_layout), alpha=0.5)
print(p)

p= ggplot(plot)+
  geom_point(aes(x=perc_error,y=ref_ratio, color=library_layout), alpha=0.5)
print(p)

p= ggplot(plot)+
  geom_point(aes(x=sum_inser,y=ref_ratio, color=library_layout), alpha=0.5)
print(p)

p= ggplot(plot)+
  geom_point(aes(x=num_genes,y=ref_ratio, color=library_layout), alpha=0.5)
print(p)

dev.off()
p= ggplot(plot)+
  geom_violin(aes(x=nError_cut,y=mode_overlap, color=library_layout))+
  geom_jitter(aes(x=nError_cut,y=mode_overlap,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="Using mean frag to calculate (mean.frag - avg_len*2)")
print(p)
p= ggplot(plot)+
  geom_violin(aes(x=nHet_cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=nHet_cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="Using mean frag to calculate (mean.frag - avg_len*2)")
print(p)

p= ggplot(plot)+
  geom_violin(aes(x=cut,y=ref_ratio, color=library_layout))+
  geom_jitter(aes(x=cut,y=ref_ratio,color=library_layout), alpha=0.5,width = 0.15)+
  labs(title="Using mean frag to calculate (mean.frag - avg_len*2)")
print(p)

dev.off()

colnames(sra_subset)

sra_subset[1:10,20:22]

