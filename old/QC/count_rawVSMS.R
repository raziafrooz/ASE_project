setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

tissue<-"Lung"
#Get the GTEx metadata:
met_gtex<-read.csv("data/GTEx_metadata.csv")
geno_met<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")
recount3_chr_mapping <- "/dcl02/lieber/ajaffe/recount-pump/recount3.alts.chromosome_mappings.tsv"
snps_gr<-readRDS("/dcs04/hansen/data/recount_genotype/biallelic_snps/biallelic_SNP_gr.rds")

k=4
  study<-met_gtex$tissue[k]
  print(k)
  load(paste0("~/hansen_lab/ASE/test_ASE/", study, "_ase_MS.rda"))
  ordered<-ase_df %>% group_by(sample_id_rep) %>%  summarize(m=median(ref_ratio)) %>% arrange(m)
  
    for(k in 1:length(unique(ordered$sample_id_rep))){
        sample_id<-ordered$sample_id_rep[k]
      
        alt_path<-met_gtex$alt[which(met_gtex$sample_id_rep==sample_id)]
        bigWig_path<-met_gtex$total[which(met_gtex$sample_id_rep==sample_id)]
        
        geno<-ase_df[ase_df$sample_id_rep==sample_id,]
        geno_gr<-makeGRangesFromDataFrame(geno, start.field = "pos", end.field = "pos", keep.extra.columns = T)
        genome(geno_gr)<-"hg38"
        ov<-findOverlaps(geno_gr,snps_gr )
        snp_gr<-snps_gr[subjectHits(ov)]
      
      
      
  #Load in bigWig file to get `coverage_count` and `filtered_snp_gr`.
  cat("Loading in: ", bigWig_path, "\n")
  
  bigwig <- tryCatch(
    {
      rtracklayer::import(bigWig_path, format = "bigwig")
    },
    error=function(cond) {
      message(paste("Error loading bigwig file: ", bigWig_path))
      message(cond)
      return(NA)
    },
    warning=function(cond) {
      message(paste("Warning loading bigwig file: ", bigWig_path))
      message(cond)
      return(NA)
    },
    finally={}
  )
  if(all(is.na(bigwig))) {
    return(NA)
  }
  overlap_loci <- GenomicRanges::findOverlaps(snp_gr, bigwig)
  bigwig_count <- bigwig$score[subjectHits(overlap_loci)]
  filter_idx <- which(bigwig_count >= 8)
  bigwig_count <- bigwig_count[filter_idx]
  filtered_snp_gr <- snp_gr[queryHits(overlap_loci)[filter_idx]]
  
  
  
  #Load in alt file to construct `alt_count`:
  #1. Decompress .zst via zstdcat in command line, save output as *.alt.temp.csv
  #2. Read in *.alt.temp.csv via vroom (fast), convert to data.table in order
  #to collapase alt counts.
  #3. Collapse to `alt_count`, and create GRanges object for intersection.
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
  ov <- GenomicRanges::findOverlaps(alt_count_gr, filtered_snp_gr)
  #Subset `alt_count_gr`, `alt_count` to the positions from SNPs, and compute
  #its own unique key.
  alt_count_gr <- alt_count_gr[queryHits(ov)]
  alt_count <- alt_count[queryHits(ov)]
  alt_key <- paste(as.character(seqnames(alt_count_gr)), start(alt_count_gr),
                   filtered_snp_gr$ref_seq[subjectHits(ov)], alt_count$alt, sep = "_")
  #Even though `alt_key` is in the same positions as `filtered_snps_key`, many of the
  #`alt_key` entries refer to alternate alleles that we are not tracking.
  #We need to match a second time, this time using the entire key.
  filtered_snps_key <- paste(as.character(seqnames(filtered_snp_gr)), start(filtered_snp_gr),
                             filtered_snp_gr$ref_seq, filtered_snp_gr$alt_seq, sep = "_")
  idx <- match(filtered_snps_key, alt_key)
  snps_key_idx <- which(!is.na(idx))
  alt_key_idx <- idx[!is.na(idx)]
  #Construct final `alt_count` relative to `filtered_snp_gr`
  final_alt_count <- rep(NA, length(filtered_snp_gr))
  final_alt_count[snps_key_idx] <- alt_count[alt_key_idx]$N
  sum(is.na(final_alt_count))
  #Compute M and S values.
  ref_count <- bigwig_count - final_alt_count
  #We sometimes have cases where the alt counts > coverage counts (bigWig):
  #the alt reads were processed to keep overlapping pair-end reads
  #whereas the coverage counts (bigWig) did not keep overlapping pair-end reads.
  #this is an ad hoc way to deal with negative counts in `ref_count`.
  
  
  geno$coverage<-bigwig_count
  geno$alt_raw<-final_alt_count
  geno$ref_raw<-ref_count
  
  saveRDS(geno, paste0("~/test/", study,sample_id,  "_raw.rda"))
    }
  
  pdf(file="~/plot/ASE/test2.pdf", width = 10, height = 4)
  for(k in 1:10){
  sample_id<-ordered$sample_id_rep[k]
  print(sample_id)
  geno<-readRDS(paste0("~/test/", study,sample_id,  "_raw.rda"))
  
 p0= ggplot(geno)+
    geom_point(aes(x=log2(ref), y=log2(alt)), alpha=0.1)+
    #geom_point(data=geno[which(geno$ref<2),],aes(x=log2(ref), y=log2(alt)), alpha=0.5, color="red")+
    geom_point(data=geno[which(geno$true_genotype==1),],aes(x=log2(ref), y=log2(alt)), alpha=0.5, color="purple")+
   geom_point(data=geno[which(geno$true_genotype==3),],aes(x=log2(ref), y=log2(alt)), alpha=0.5, color="blue")+
   labs(caption=paste0("old= ", median(geno$ref_ratio), ", new=",
                       median(geno$ref_ratio[-which(geno$true_genotype!=2)]) ),
        title=paste0("10 samples in Lung (ordered low to high ref_ratio)",sample_id) ,
        subtitle= "purple is true genotype = homo_ref, blue is homo_alt" )
 print(p0)
  }
  dev.off()

  
  pdf(file="~/plot/ASE/test3.pdf", width = 10, height = 4)
  for(k in 1:10){
    sample_id<-ordered$sample_id_rep[k]
    print(sample_id)
    geno<-readRDS(paste0("~/test/", study,sample_id,  "_raw.rda"))
    
    p0= ggplot(geno)+
      geom_point(aes(x=log2(ref), y=log2(alt)), alpha=0.1)+
      #geom_point(data=geno[which(geno$ref<2),],aes(x=log2(ref), y=log2(alt)), alpha=0.5, color="red")+
      geom_point(data=geno[which(geno$ref_ratio<0.25),],aes(x=log2(ref), y=log2(alt)), alpha=0.5, color="purple")+
      geom_point(data=geno[which(geno$ref_ratio>0.75),],aes(x=log2(ref), y=log2(alt)), alpha=0.5, color="blue")+
      labs(caption=paste0("old= ", median(geno$ref_ratio), ", new=",
                          median(geno$ref_ratio[-which(geno$ref_ratio<0.2 | geno$ref_ratio> 0.7)]) ),
           title=paste0("10 samples in Lung (ordered low to high ref_ratio)",sample_id) ,
           subtitle= "purple is ref_ratio<0.25, blue is ref_ratio>0.75" )
    print(p0)
  }
  dev.off()
  
  all.equal(geno_10$ref, geno_10$ref_raw)
  k=4
  for(k in 1:10){
    sample_id<-ordered$sample_id_rep[k]
    print(sample_id)
    geno<-readRDS(paste0("~/test/", study,sample_id,  "_raw.rda"))
    geno_gr<-makeGRangesFromDataFrame(geno, start.field = "pos", end.field = "pos", keep.extra.columns = T)
    
    ov<-findOverlaps(black_list,geno_gr)
    print(median(geno_gr$ref_ratio))
    print(median(geno_gr$ref_ratio[-subjectHits(ov)]))
  }
  
  black_list
  bed_graph
  ov<-findOverlaps(black_list,geno_gr)
  geno_gr[subjectHits(ov)]
  
  geno_gr[]
  
  pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
  for(k in 1:10){
    sample_id<-ordered$sample_id_rep[k]
    print(sample_id)
    geno<-readRDS(paste0("~/test/", study,sample_id,  "_raw.rda"))
    
  p0=ggplot(geno)+
    geom_histogram(aes(ref_ratio), fill="salmon",alpha=0.5)+
    geom_histogram(data=geno[geno$true_genotype!=2,], aes(ref_ratio), fill="purple", alpha=0.5)+
    labs(title=sample_id,
         subtitle="salmon = all snps, purple= genotyping errors (!=2)")
  print(p0)
  }
dev.off()  
  
geno<-geno %>% mutate(cut=cut_number(total,5))

pdf(file="~/plot/ASE/test1.pdf", width = 10, height = 4)
  p0=ggplot(geno)+
    geom_histogram(aes(ref_ratio,fill=cut),alpha=0.5)+
    facet_wrap(vars(cut))
  
  print(p0)

dev.off()  
  #[1] 0.4897593
  mean(geno$ref_ratio[geno$total>15 & geno$true_genotype==2] )
  mean(geno$ref_ratio[geno$ref_ratio>0.12 & geno$ref_ratio<0.87])
  
  
  geno[geno$ref<2,c("chr","pos","total","ref","alt")]
  geno_gr<-makeGRangesFromDataFrame(geno, start.field = "pos", end.field = "pos", keep.extra.columns = T)
  
  
  
  
  #--------------------------------------------------------
  #Variant annotation
  #--------------------------------------------------------
  # if (!require("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # 
  # BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
  library(VariantAnnotation)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  #seqlevels(vcf) <- "chr22"
  #rd <- rowRanges(vcf)
  loc <- locateVariants(geno_gr, txdb, AllVariants())
  head(loc, 6)
  ov<-findOverlaps(geno_gr,loc)
  ov<-ov[!duplicated(queryHits(ov))]
  geno$annotation<-NA
  geno$annotation[queryHits(ov)]<-as.character(loc$LOCATION[subjectHits(ov)])
  which(is.na(geno$annotation))
  loc[which(loc$QUERYID ==2)]
  geno[2117,]
  
  geno %>% group_by(annotation) %>% summarize(median(ref_ratio), sum(true_genotype!=2), n())
  
  
  
  pdf(file="~/plot/ASE/test3.pdf", width = 10, height = 4)
    p0= ggplot(geno %>% filter(annotation%in% c("intron", "NA","intergenic" )))+
      geom_point(aes(x=log2(ref), y=log2(alt), color=annotation), alpha=0.5)
  print(p0)
  p0= ggplot(geno)+
    geom_point(aes(x=log2(ref), y=log2(alt)), alpha=0.5)
  print(p0)
  dev.off()
  
  geno %>% filter(!annotation %in% c("intron", "NA","intergenic" )) %>% summarize(median(ref_ratio))
  
  median(geno$ref_ratio)
  
  
  #--------------------------------------------------------
  #Alexis Battle SRA
  #--------------------------------------------------------
  
  sra_battle<-read.delim("data/battle.tsv", header=T, sep="\t")
  recount3_metadata<-fread("/dcs04/hansen/data/recount_genotype/PCA/SRA/Recount3_metadata.tsv", header= T, sep = "\t",quote="")
  recount3_metadata<-recount3_metadata[,1:6]
  sra<-readRDS("data/use_allSRA_QC.rds")
  ?read.table
  sra$exp_acc<-recount3_metadata$experiment_acc[match(sra$sample_id,recount3_metadata$external_id)]
  sra<- sra %>% relocate(exp_acc, .after=sample_id)
  sra[35:40,]
  recount3_metadata[1:4,]
  recount3_metadata[1,]
  sra_battle[1:5,]
  read.tabl
  plot<-sra[which(sra$study %in% sra_battle$study),]
  sra$tissue<-sra_battle$tissue[match(sra$exp_acc, sra_battle$sample)]
  sra$tissue.category<-sra_battle$tissue.category[match(sra$exp_acc, sra_battle$sample)]
  sra$disease.category<-sra_battle$disease.category[match(sra$exp_acc, sra_battle$sample)]
  sra$cancer<-sra_battle$cancer[match(sra$exp_acc, sra_battle$sample)]
  
  sra_battle[1,]
  sra_battle[which(sra_battle$sample=="SRX683800"),]
  sra_filt<-sra[which(sra$exp_acc %in% sra_battle$sample),]
  duplicated(sra_filt$exp_acc)[1:100]
  sra_filt[95:100,]
  
  xx<-sra_filt[sra_filt$exp_acc=="SRX683800",]
  
  load(paste0("~/hansen_lab/ASE/test_ASE/SRR1554542_ase.rda")  )
  
  run2 <- ase_df
  rm(ase_df)
  all_runs<-run1 %>% bind_rows(run2) %>%  group_by(chr, pos) %>%  summarize(total=sum(total), ref=sum(ref), alt=sum(alt))
  all_runs$ref_ratio<-all_runs$ref/all_runs$total
median(all_runs$ref_ratio[which(all_runs$total>20)])  

unique(sra$tissue.category)[1:40]


library(ggridges)

pdf(file="~/plot/ASE/sra_alexis.pdf", width = 10, height = 4)

ggplot(plot,aes(ref_ratio))+
  geom_histogram(bins = 100)+ 
  geom_vline(xintercept = 0.5, color="red")+ geom_vline(xintercept = 0.55, color="purple")+ geom_vline(xintercept = 0.45, color="blue")+
  labs(title="SRA studies that are selected by Alexis Battle (~90K samples)", 
       subtitle="blue=0.45, red=0.5, purple=0.55")+
  xlim(c(0.2,0.8))
dev.off()


  
  library(recount3)
  
  
  human_projects <- available_projects()
  proj_info <- subset(
    human_projects,
    project =="GTEX-13OW7-0226-SM-5MR3N.1"  & project_type == "data_sources"
  )
  rse <- create_rse(proj_info)
  ## Create a RangedSummarizedExperiment (RSE) object at the gene level
  #rse_gene_SRP009615 <- create_rse(proj_info)
  #colData(rse_gene_SRP009615)[1,30:40]
  #quantile(rse_gene_SRP009615$sra.run_total_bases)
  url<-locate_url(
    "ADIPOSE_SUBCUTANEOUS",
    "data_sources/gtex",
    type = "metadata")
  
  x <-utils::read.delim(file_retrieve(url[1], verbose = FALSE))
  x[1,89]
  x[1,1:5]
  x[x$external_id==sample_id,]
  
  unique(x$external_id)[20:30]
  
  
  
  
  #----------------
  ov <- findOverlaps(filtered_snp_gr,alt_count_gr)
  which(duplicated(queryHits(ov)))[76]
  ov[560:563]
  geno[170,c("chr", "pos", "total", "ref","alt","noise")]
  alt_count[c(2509209,2509215),]
  which(geno$ref<=)
  which(which(duplicated(queryHits(ov))) %in% which(queryHits(ov) %in% which(geno$ref<=4)))
  
  noise<-alt_count[!(alt_key %in% filtered_snps_key),]
  noise_gr <- GenomicRanges::GRanges(seqnames = noise$chr,
                                     ranges = IRanges(noise$pos,
                                                      noise$pos))
  ov <- GenomicRanges::findOverlaps(noise_gr, filtered_snp_gr)
  geno$noise<- NA
  geno$noise[subjectHits(ov)]<-noise$N[queryHits(ov)]
  #idx <- match(filtered_snps_key, alt_key)
  #snps_key_idx <- which(!is.na(idx))
  #alt_key_idx <- idx[!is.na(idx)]
  #Construct final `alt_count` relative to `filtered_snp_gr`
  #final_alt_count <- rep(NA, length(filtered_snp_gr))
  #final_alt_count[snps_key_idx] <- alt_count[alt_key_idx]$N
  geno<-geno %>% rowwise %>% mutate(min_allel=min(ref,alt))
  geno$noise_perc<-geno$noise/geno$total*100
  median(geno$ref_ratio[-which(geno$noise_perc>=5)])
  pdf(file="~/plot/ASE/test5.pdf", width = 10, height = 4)
  ggplot(geno)+
    geom_point(aes(x=log2(ref), y=log2(alt)), alpha=0.1)+
    geom_point(data=geno[which(geno$noise_perc>=5),],aes(x=log2(ref), y=log2(alt)), alpha=0.5, color="purple")
  ggplot(geno)+
    geom_point(aes(x=total, y=min_allel), alpha=0.1)+
    geom_point(data=geno[which(geno$noise_perc>=5),], aes(x=total, y=min_allel), color="purple")+
    xlim(c(0,500))+
    ylim(c(0,100))
  dev.off()
  
  geno4<-geno
  df<-data.frame(ref_ratio=rep(NA,20), total_noise=rep(NA,20), noise_alot=rep(NA,20))
  df$ref_ratio[6]<-median(geno6$ref_ratio)
  df$total_noise[6]<-sum(!is.na(geno6$noise))
  df$noise_alot[6]<-sum(geno6$noise_perc>=5)
  #----------------
  
  
  for(k in 7:20){
    sample_id<-ordered$sample_id_rep[k]
    
    alt_path<-met_gtex$alt[which(met_gtex$sample_id_rep==sample_id)]
    bigWig_path<-met_gtex$total[which(met_gtex$sample_id_rep==sample_id)]
    
    geno<-ase_df[ase_df$sample_id_rep==sample_id,]
    geno_gr<-makeGRangesFromDataFrame(geno, start.field = "pos", end.field = "pos", keep.extra.columns = T)
    genome(geno_gr)<-"hg38"
    ov<-findOverlaps(geno_gr,snps_gr )
    snp_gr<-snps_gr[subjectHits(ov)]
    
    
    
    #Load in bigWig file to get `coverage_count` and `filtered_snp_gr`.
    cat("Loading in: ", bigWig_path, "\n")
    
    bigwig <- tryCatch(
      {
        rtracklayer::import(bigWig_path, format = "bigwig")
      },
      error=function(cond) {
        message(paste("Error loading bigwig file: ", bigWig_path))
        message(cond)
        return(NA)
      },
      warning=function(cond) {
        message(paste("Warning loading bigwig file: ", bigWig_path))
        message(cond)
        return(NA)
      },
      finally={}
    )
    if(all(is.na(bigwig))) {
      return(NA)
    }
    overlap_loci <- GenomicRanges::findOverlaps(snp_gr, bigwig)
    bigwig_count <- bigwig$score[subjectHits(overlap_loci)]
    filter_idx <- which(bigwig_count >= 8)
    bigwig_count <- bigwig_count[filter_idx]
    filtered_snp_gr <- snp_gr[queryHits(overlap_loci)[filter_idx]]
    
    
    
    #Load in alt file to construct `alt_count`:
    #1. Decompress .zst via zstdcat in command line, save output as *.alt.temp.csv
    #2. Read in *.alt.temp.csv via vroom (fast), convert to data.table in order
    #to collapase alt counts.
    #3. Collapse to `alt_count`, and create GRanges object for intersection.
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
    ov <- GenomicRanges::findOverlaps(alt_count_gr, filtered_snp_gr)
    #Subset `alt_count_gr`, `alt_count` to the positions from SNPs, and compute
    #its own unique key.
    alt_count_gr <- alt_count_gr[queryHits(ov)]
    alt_count <- alt_count[queryHits(ov)]
    alt_key <- paste(as.character(seqnames(alt_count_gr)), start(alt_count_gr),
                     filtered_snp_gr$ref_seq[subjectHits(ov)], alt_count$alt, sep = "_")
    #Even though `alt_key` is in the same positions as `filtered_snps_key`, many of the
    #`alt_key` entries refer to alternate alleles that we are not tracking.
    #We need to match a second time, this time using the entire key.
    filtered_snps_key <- paste(as.character(seqnames(filtered_snp_gr)), start(filtered_snp_gr),
                               filtered_snp_gr$ref_seq, filtered_snp_gr$alt_seq, sep = "_")
    
    
    noise<-alt_count[!(alt_key %in% filtered_snps_key),]
    noise_gr <- GenomicRanges::GRanges(seqnames = noise$chr,
                                       ranges = IRanges(noise$pos,
                                                        noise$pos))
    ov <- GenomicRanges::findOverlaps(noise_gr, filtered_snp_gr)
    geno$noise<- NA
    geno$noise[subjectHits(ov)]<-noise$N[queryHits(ov)]
    
    df$ref_ratio[k]<-median(geno$ref_ratio)
    df$total_noise[k]<-sum(!is.na(geno$noise))
    df$noise_alot[k]<-sum(geno$noise_perc>=5)
  }
  
  df<-df %>% mutate(cut=cut_number(total_noise,5))
  pdf(file="~/plot/ASE/test1.pdf", width = 10, height = 4)
  p0=ggplot(df)+
    geom_histogram(aes(ref_ratio,fill=cut),alpha=0.5)
  
  print(p0)
  
  dev.off()  
  
  #---------------------------------------------------------------
  
  plot<-data.frame(g_error8=rep(NA, 30), g_error10=rep(NA, 30),g_error15=rep(NA, 30),g_error20=rep(NA, 30),g_error30=rep(NA, 30),
                   total_8=rep(NA, 30), total_10=rep(NA, 30), total_15=rep(NA, 30), total_20=rep(NA, 30),total_30=rep(NA, 30))
  for(k in 1:30){
    sample_id<-ordered$sample_id_rep[k]
    print(sample_id)
    geno<-readRDS(paste0("~/test/", study,sample_id,  "_raw.rda"))
    
  for(i in 1:5){
    lim<-c(8,10,15,20,30)
    xx<-geno %>% filter(total>=lim[i])
    plot[k,i]<-sum(xx$true_genotype!=2)/nrow(xx)
    plot[k,i+5]<-median(xx$ref_ratio)
  }
  }
  k=20
  i=1
  plot$sample_id<-ordered$sample_id_rep[1:30]
  plot0<- plot %>% pivot_longer(!sample_id, names_to="cut_off", values_to="ref_ratio")
  pdf(file="~/plot/ASE/test.pdf", width = 10, height = 4)
  p0= ggplot(plot0)+
    geom_point(aes(x=sample_id, y=ref_ratio, color=cut_off), alpha=0.8)
  print(p0)
  dev.off()
  