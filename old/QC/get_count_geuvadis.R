# In this analysis I am using the ASE data from GTEx
#This data was downloaded by Nick during rotations from dbGap: V8 not WASP corrected :location /dcl01/hansen/data/gtex_ase
#The GTEx paper recommends using the ase WASP corrected (allelic mapping error). Afrooz downloaded this on 08/09/2023 and is here /dcl01/hansen/data/arazi/ASE/dbGap

###For now start the analysis with ase not WASP corrected:

setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)


recount3_chr_mapping <- "/dcl02/lieber/ajaffe/recount-pump/recount3.alts.chromosome_mappings.tsv"
snp_gr<-readRDS("/dcs04/hansen/data/recount_genotype/biallelic_snps/biallelic_SNP_gr.rds")

#met<-read.csv("data/GTEx_metadata.csv")
#geno_met<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")

geu<-read.csv("~/ASE/data/Geuvadis_metadata.csv")
geu_geno_meta<- "/dcs04/hansen/data/recount_genotype/pipeline/ERP001942/predict_genotype_accuracy/"

#########################
#start with kidney medulla since it only has 4 samples
#########################
#tissue<-c("Liver", "Lung", "Stomach", "Pancreas")
#tissue_ab<-c("LIVER", "LUNG","STMACH","PNCREAS")

for(k in 1:length(geu$sample_id_rep)){
  #tissue_name<-tissue[k]
  sample_id<-geu$sample_id_rep[k]
  study<-geu$study[k]
  print(k)
  
  if(!file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, ".rds")))
  {
    geno_file_path<-paste0(geu_geno_meta, sample_id,"_predGenotypes_w_accuracy.csv.gz" )
    alt_path<-geu$alt[which(geu$sample_id_rep==sample_id)]
    bw_path<-geu$total[which(geu$sample_id_rep==sample_id)]
    print(sample_id)
    geno_data<-as_tibble(read.csv(geno_file_path))
    #geno_data<-readRDS(geno_met$allGenotypesOutput[geno_met$study==tissue_name])
    #for(i in 1:length(unique(geno_data$sample_id_rep))){
     # print(i)
      #sam<-unique(geno_data$sample_id_rep)[i]
      ##start the analysis with 1 indv, only select hetrozyogus locations, total read cound > 8)
      geno<-geno_data %>% filter(pred_genotype==2, coverage>=8) 
      geno_gr<-makeGRangesFromDataFrame(geno, start.field = "start", end.field = "start", keep.extra.columns = T)
      genome(geno_gr)<-"hg38"
      
      
      #############
      #Read in the alt file
      ##############
      temp="~/test"
      temp_altFile <- paste0(temp, sample_id, ".alt.temp.csv")
      system(paste0("zstdcat ", alt_path, " > ", temp_altFile))
      alt <- fread(temp_altFile, select = c(1, 2, 4))
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
      alt_count_gr <- GRanges(seqnames = alt_count$chr,
                              ranges = IRanges(alt_count$pos,
                                               alt_count$pos))
      
      ov <- findOverlaps(alt_count_gr,snp_gr)
      
      #################################################
      #this is not exclusively biallelic!!!!!!!!!!!! come back to it and filter
      #################################################
      #Subset `alt_count_gr`, `alt_count` to the positions from SNPs, and compute 
      #its own unique key.
      alt_count_gr <- alt_count_gr[queryHits(ov)]
      alt_count <- alt_count[queryHits(ov)]
      alt_key <- paste(as.character(seqnames(alt_count_gr)), start(alt_count_gr), 
                       snp_gr$ref_seq[subjectHits(ov)], alt_count$alt, sep = "_")
      #Even though `alt_key` is in the same positions as `filtered_snps_key`, many of the 
      #`alt_key` entries refer to alternate alleles that we are not tracking. 
      #We need to match a second time, this time using the entire key.
      ov <- findOverlaps(geno_gr,snp_gr)
      filtered_snps_key <- paste(as.character(seqnames(geno_gr)), start(geno_gr), 
                                 snp_gr$ref_seq[subjectHits(ov)], snp_gr$alt_seq[subjectHits(ov)], sep = "_")
      idx <- match(filtered_snps_key, alt_key)
      snps_key_idx <- which(!is.na(idx))
      alt_key_idx <- idx[!is.na(idx)]
      #Construct final `alt_count` relative to `filtered_snps_gr`
      final_alt_count <- rep(0, length(geno_gr))
      final_alt_count[snps_key_idx] <- alt_count[alt_key_idx]$N
      
      
      geno_gr$alt_count<-final_alt_count
      rm(final_alt_count)
      #We sometimes have cases where the alt counts > coverage counts (bigWig): 
      #the alt reads were processed to keep overlapping pair-end reads
      #whereas the coverage counts (bigWig) did not keep overlapping pair-end reads. 
      #this is an ad hoc way to deal with negative counts in `ref_mtx`.
      n_dupli<-length(geno_gr[geno_gr$alt_count > geno_gr$coverage])
      print(paste0("duplicated_counts = ", n_dupli))
      geno_gr[geno_gr$alt_count > geno_gr$coverage] <-NULL ###for now just exclude them
      
      #Get ref count:
      geno_gr$ref_count<-geno_gr$coverage-geno_gr$alt_count
      
      
      ase_df<-tibble(chr=as.character(seqnames(geno_gr)),pos=start(geno_gr),
                      sample_id=sample_id,
                      ref=geno_gr$ref_count, 
                      alt= geno_gr$alt_count, 
                      total=geno_gr$coverage,
                      pred_genotype=geno_gr$pred_genotype,
                      num_dupl=n_dupli)
      
      IS.MEDIAN<-median(ase_df$ref/ase_df$total)
      print(paste0("Median is ", IS.MEDIAN ))
      warning(if(!0.4<IS.MEDIAN &IS.MEDIAN <0.6){print("Median is not 0.5")}) 
      
      
      
    saveRDS(ase_df,file = paste0("~/hansen_lab/ASE/test_ASE/", sample_id, ".rds"))
  }
}


library(qvalue)

#Gtex paper suggests removing vriants in HLA genes:
#his BED file contains genomic positions that we have identified as either showing bias in simulations or having a UCSC mappability score < 50. 
#Variants that fall into these positions are used for phasing, but not for generating haplotypic counts to avoid problems with mapping bias.
#https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
#Downloaded blacklist of snps from gtex (https://github.com/secastel/phaser/blob/master/phaser/README.md)
bad_snp<-read.table("~/plot/ASE/hg38_haplo_count_blacklist.chr.bed", sep="\t")
bad_snp_gr<-makeGRangesFromDataFrame(bad_snp,seqnames="V1",start.field ="V2",end.field = "V3")


k=1
for(k in 1:length(geu$sample_id_rep)){
  #tissue_name<-tissue[k]
  sample_id<-geu$sample_id_rep[k]
  study<-geu$study[k]
  print(k)
  if(!file.exists(file=paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda") ))
  {
    print(sample_id)
    ase_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, ".rds"))
    ase_df$ref_ratio<- ase_df$ref/ase_df$total
      ase_filt_gr<-makeGRangesFromDataFrame(ase_df,seqnames="chr",start.field ="pos",end.field = "pos", keep.extra.columns = T)
      
      #Remove blacklist snps from our granges:
      ov<-findOverlaps(ase_filt_gr,bad_snp_gr)
      ase_filt_gr<-ase_filt_gr[-unique(queryHits(ov))]
      print(paste0("number of blacklist: ", length(unique(queryHits(ov)))))
      
      #obtain the median of the ref ratio to be used as the null p value
      is.median<-median(ase_filt_gr$ref_ratio)
      print(paste0("Median is: ", is.median))
      
      ase_filt_gr$p_val = apply(as.data.frame(ase_filt_gr)[,c("ref","alt")], 1, function(x) binom.test(x[1],(x[1]+x[2]),p=is.median)$p.value)
      # perform multiple testing correction with FDR
      ase_filt_gr$q_val = p.adjust(ase_filt_gr$p_val, method = "fdr")
      
     
        ase_all<- as_tibble(ase_filt_gr)
      
    #order the df and save the data:
    colnames(ase_all)[1:2]<-c("chr","pos")
    ase_all[,3:5]<-NULL
    save(ase_all, file=paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda") )
  }
}

#----------------------------------------
# Try without blacklist:
#----------------------------------------

for(k in 1:length(geu$sample_id_rep)){
  #tissue_name<-tissue[k]
  sample_id<-geu$sample_id_rep[k]
  study<-geu$study[k]
  print(k)
    print(sample_id)
    ase_df<-readRDS(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, ".rds"))
    ase_df$ref_ratio<- ase_df$ref/ase_df$total
    #obtain the median of the ref ratio to be used as the null p value
    is.median<-median(ase_df$ref_ratio)
    print(paste0("Median is: ", is.median))
    
    #ase_df$p_val = apply(ase_df[,c("ref","alt")], 1, function(x) binom.test(x[1],(x[1]+x[2]),p=is.median)$p.value)
    # perform multiple testing correction with FDR
    #ase_df$q_val = p.adjust(ase_df$p_val, method = "fdr")
    
    
    ase_all<- as_tibble(ase_df)
    if(k==1){
      ase_all_geu_noBlack<-ase_all
    }
    ase_all_geu_noBlack<-rbind(ase_all_geu_noBlack,ase_all)
}




#-------------------
#plot
#-------------------

for(k in 1:length(geu$sample_id_rep)){
  #tissue_name<-tissue[k]
  sample_id<-geu$sample_id_rep[k]
  study<-geu$study[k]
  print(k)
  if(file.exists(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda")))
  {
    load(paste0("~/hansen_lab/ASE/test_ASE/", sample_id, "_ase.rda"))
   
    if(k==1){
      ase_all_geu<-ase_all
    }
    ase_all_geu<-rbind(ase_all_geu,ase_all)
  }
}

pdf(file="~/plot/ASE/geuvadis_rawCounts.pdf", width = 10, height = 4)

ggplot(data=ase_all_geu, aes(x=sample_id,y=mean(ref_ratio)))+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")

dev.off()
df<-ase_all_geu %>% group_by(sample_id) %>% summarize(M=mean(ref_ratio))
df<-ase_all_geu_noBlack %>% group_by(sample_id) %>% summarize(M=mean(ref_ratio))
pdf(file="~/plot/ASE/geuvadis_rawCounts_noBlacklist.pdf", width = 10, height = 4)

ggplot(data=df, aes(x=1:462,y=M))+
  geom_point()+#outlier.shape = NA)+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  labs(x="sample ID (1:462)", 
       title= "All Geuvadis samples: calculated mean ref_ratio for each sample",
       subtitle="Each dot is a sample's mean ref_ratio")

ggplot(data=ase_all_geu %>% filter(sample_id %in% unique(ase_all_geu$sample_id)[1:200]), aes(x=sample_id,y=ref_ratio))+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  labs(x="sample ID (1:200)",
       title= "First half of Geuvadis samples: calculated mean ref_ratio for each sample",
       subtitle="Each box is a sample")

ggplot(data=ase_all_geu%>% filter(sample_id %in% unique(ase_all_geu$sample_id)[201:462]), aes(x=sample_id,y=ref_ratio))+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
  labs(x="sample ID (201:462)",
       title= "Second half of Geuvadis samples: calculated mean ref_ratio for each sample",
       subtitle="Each box is a sample")

dev.off()





pdf(file="~/plot/ASE/mono_allelic_filter_geu.pdf", width = 10, height = 4)
for(i in unique(ase_all_geu$sample_id)[1:10]){
  print(i)
  x<-ase_all_geu %>% filter(sample_id==i)%>% rowwise() %>% mutate(mine_allele=min(alt,ref))
  p=ggplot(x)+
    geom_point(aes(x=total, y=mine_allele))+xlim(0,100)+ylim(0,100)+
    labs(title="First 10 Geuvadis samples: mono-allele check")
  print(p)
}
dev.off()


pdf(file="~/plot/ASE/geuvadis_coverage.pdf", width = 10, height = 4)
for(i in unique(ase_all_geu$sample_id)[1:10]){
  x<-ase_all_geu %>% filter(sample_id==i)
  p=ggplot(x)+
    geom_histogram(aes(total))+
    xlim(c(0,1000))+
    labs(title="Geuvadis coverage distribution",
         subtitle=paste0("sample=",unique(x$sample_id)[1]))
  print(p)
}
dev.off()

#--------------
#Plot population
#--------------
#Download population for geuvadis from here: https://www.internationalgenome.org/data-portal/sample

pop<-read.csv("data/geuvadis_pop.tsv", sep = "\t")
ase_all_geu$indv_id<-geu$individual_id[match(ase_all_geu$sample_id,geu$sample_id_rep)]
ase_all_geu$pop<-pop$Population.code[match(ase_all_geu$indv_id,pop$Sample.name)]
ase_all_geu$super_pop<-pop$Superpopulation.code[match(ase_all_geu$indv_id,pop$Sample.name)]

plot_df<-ase_all_geu %>% group_by(sample_id) %>% summarize(ref_ratio= median(ref_ratio))

plot_df$indv_id<-geu$individual_id[match(plot_df$sample_id,geu$sample_id_rep)]
plot_df$pop<-pop$Population.code[match(plot_df$indv_id,pop$Sample.name)]
plot_df$super_pop<-pop$Superpopulation.code[match(plot_df$indv_id,pop$Sample.name)]

pdf(file="~/plot/ASE/geuvadis_pop.pdf", width = 10, height = 4)
for (i in unique(ase_all_geu$pop)) {
  p=ggplot(ase_all_geu %>% filter(pop==i))+
    geom_boxplot(aes(x=sample_id, y=ref_ratio),outlier.shape = NA)+
    geom_hline(yintercept = 0.5, color="red",linetype="dotdash")+
    labs(title=paste0("Geuvadis sample from population = ",i))
  print(p)
}

ggplot(plot_df)+
  geom_point(aes(x=pop, y=ref_ratio, color=super_pop))+
  geom_hline(yintercept = 0.5, color="red",linetype="dotdash")
dev.off()

pdf(file="~/plot/ASE/geuvadis_test.pdf", width = 10, height = 4)
d=ase_all_geu %>% filter(sample_id== "ERR188482")
ggplot(d)+
  geom_point(aes(x=log2(ref), y=log2(alt)))+
  geom_point(data = d %>% filter(q_val<0.05), aes(x=log2(ref), y=log2(alt)), color="red")

dev.off()

d %>% filter(q_val<0.05, log2(alt)<5,log2(ref)>7)

