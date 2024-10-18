###################
#Geuvadis metadata:
###################

#change this metadata to new locations: https://github.com/raziafrooz/recount_genotype/blob/main/study_metadata/geuvadis/geuvadis_metadata.csv
geu<-read.csv("~/test/geuvadis_metadata.csv")

#define constant varaibles of where bigwigs and alts are located. 
file_path <- "/dcs04/hansen/data/recount3/geuvadis"
setwd("/dcs04/hansen/data/recount3/geuvadis/")
#get bigwig files
bigwig_files <- system("find . -name *.all.bw", intern = T)
bigwig_files <-sub('.', '', bigwig_files)

path_total_bw <- paste0(file_path, bigwig_files)

##get alt file:
alt_files <- system("find . -name *bamcount_nonref.csv.zst", intern = T)
alt_files <-sub('.', '', alt_files)

path_alt <- paste0(file_path, alt_files)

sample <- sapply(strsplit(bigwig_files, "[/]"), function(x) x[3])
df<-data.frame(sample, path_total_bw,path_alt)

geu$total<-df$path_total_bw[match(geu$sample_id_rep,df$sample)]
geu$alt<-df$path_alt[match(geu$sample_id_rep,df$sample)]
write.csv(geu,"~/ASE/data/Geuvadis_metadata.csv")
