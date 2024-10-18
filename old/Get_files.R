setwd("~/ASE/")
library(tidyverse)

#define constant varaibles of where bigwigs and alts are located. 
###This should be the final path but as of 09/05/2023 the files are not here yet
#bigwig_path <- "/dcs04/hansen/data/recount3/gtex_monorail_output/bigwigs/" 
#alt_path <- "/dcs04/hansen/data/recount3/gtex_monorail_output/alts/" 
###They are temporarly here
bigwig_path <- "/dcl02/lieber/ajaffe/trash_jhpce_org/recount-pump/gtex_monorail_output/bigwigs/"
alt_path <- "/dcl02/lieber/ajaffe/trash_jhpce_org/recount-pump/gtex_monorail_output/alts/"
tissue_info_path <- "/users/swang1/gtex_count/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

#get bigwig files
bigwig_files <- system(paste0("/bin/ls ", bigwig_path, " | grep .all.bw"), intern = T)
path_total_bw <- paste0(bigwig_path, bigwig_files)
#get alt files
sample_rep <- sapply(strsplit(bigwig_files, "[.]"), function(x) x[2])
sample <- sub(".all.bw", "", bigwig_files)
sample_id <- sapply(strsplit(bigwig_files, "[.]"), function(x) x[1])
path_alts_bw <- paste0(alt_path, sample, ".bamcount_nonref.csv.zst")

#make a data table with all the sample names and file locations 
metadata <- data.frame(sample_id = sample_id, rep_id=sample_rep, 
                       total= path_total_bw, alt = path_alts_bw)

#get tissue names
sample_attribute <- read.delim(tissue_info_path, sep = "\t", 
                               stringsAsFactors=F)
metadata$tissue <- as.character(sapply(1:nrow(metadata), 
                                       function(x) sample_attribute$SMTSD[sample_attribute$SAMPID == metadata$sample_id[x]]))

#get individual ID
metadata$individual_id <- unlist(lapply(strsplit(metadata$sample_id, "-"), 
                                        function(x) paste0(x[1], "-", x[2])))

#add column for combination of sample_id and replicates:
metadata$sample_id_rep <- paste0(metadata$sample_id, ".", 
                                 metadata$rep_id)

#filter out missing tissue samples
metadata <- metadata[metadata$tissue != "character(0)" ,]

#change tissue name to have no "-" or " " or "(" or ")" characters.
#those characters will crash Snakemake. 
metadata$tissue <- gsub("\\s*\\([^\\)]+\\)", "", 
                        metadata$tissue) #removes ()
metadata$tissue <- trimws(metadata$tissue)
metadata$tissue <- gsub(" - ", "_", metadata$tissue)
metadata$tissue <- gsub(" ", "_", metadata$tissue)

#create study column: the column Snakemake will group each 
#study/tissue in its data processing.
metadata$study <- metadata$tissue

#now, generate metadatas for snakemake
write.csv(metadata, "data/GTEx_metadata.csv", quote = F, row.names = F)

