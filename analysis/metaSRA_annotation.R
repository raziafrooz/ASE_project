library(RSQLite)
library(DBI)
library(tidyverse)
library(DO.db)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(data.table)
setwd("~/ASE/")

#------------------------------
#Get Cancer ID from Disease Ontology
#------------------------------
xx <- as.list(DOOFFSPRING)
doid_tb<-tibble(doid=names(xx))
#Cancer terms: found maually from here https://disease-ontology.org/do
#sarcoma parent DOID:1115,disease of cellular proliferation "DOID:14566", connective tissue DOID:201 and DOID:65,
#carcenoma DOID:305, cell type cancer DOID:0050687, organ system cancer DOID:0050686,
#hematopoietic system disease DOID:74


yy<-unique(c(xx[[which(doid_tb$doid=="DOID:162")]],
      xx[[which(doid_tb$doid=="DOID:1115")]],
      xx[[which(doid_tb$doid=="DOID:1749")]],
      xx[[which(doid_tb$doid=="DOID:14566")]],
      xx[[which(doid_tb$doid=="DOID:201")]],
      xx[[which(doid_tb$doid=="DOID:65")]],
      xx[[which(doid_tb$doid=="DOID:305")]],
      xx[[which(doid_tb$doid=="DOID:0050687")]],
      xx[[which(doid_tb$doid=="DOID:0050686")]],
      xx[[which(doid_tb$doid=="DOID:162")]]))

# FirstTenDOBimap <- DOTERM
# xx <- as.list(FirstTenDOBimap)
# 
# 
# ds_term<-unlist(lapply(xx,function(y) Term(y)))
# doid_tb<-tibble(doid=names(xx), temrs=ds_term)
# 
# toMatch<-c("cancer", "tumor","carcinoma", "sarcomas","leukemia","oma")
# index <- grep(paste(toMatch,collapse="|"),doid_tb$temrs, ignore.case = TRUE)
# doid_tb$cancer<-"none_cancer"
# doid_tb$cancer[index]<-"cancer"
# yy<-unique(doid_tb$doid[which(doid_tb$cancer=="cancer")])



#------------------------------
#Get sample annotation from metaSRA
#------------------------------

#download the sqlite latest version from https://metasra.biostat.wisc.edu/supportpages/download.html
#load it in:
connection<-dbConnect(drv=RSQLite::SQLite(), dbname="/dcs07/hansen/data/recount_ASE/data/metasra.v1-8.sqlite")
#dbListTables(connection)
cell_type<-tbl(connection, "sample_type") 

rs_cell <- dbSendQuery(connection, "SELECT * FROM sample_type")
cell_type <- dbFetch(rs_cell)

rs <- dbSendQuery(connection, "SELECT * FROM mapped_ontology_terms")
rs1 <- dbFetch(rs)
cancer_samp<- rs1 %>% 
  filter(term_id %in% yy) %>%
  dplyr::select(sample_accession) %>% unique()

saveRDS(rs1, file="/dcs07/hansen/data/recount_ASE/data/sra_ontology_term.rds")
saveRDS(cancer_samp, file="/dcs07/hansen/data/recount_ASE/data/cancer_annot.rds")
saveRDS(cell_type, file="/dcs07/hansen/data/recount_ASE/data/cell_type_annot.rds")

#------------------------------------------
#MetaSRA is missing annotation for some samples:
#do a word search:
#------------------------------------------
ase_df<-fread("/dcs07/hansen/data/recount_ASE/data/sra_ASE_stat.csv.gz")
recount3_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/metadata/Recount3_metadata.tsv", header= T, sep = "\t",quote="")


ase_df$sample_attributes<-recount3_metadata$sample_attributes[match(ase_df$experiment_acc,recount3_metadata$experiment_acc)]
ase_df$sample_acc<-recount3_metadata$sample_acc[match(ase_df$experiment_acc,recount3_metadata$experiment_acc)]

not_found<-ase_df[(is.na(match(ase_df$sample_acc,rs1$sample_accession))),]


toMatch<-c("cancer", "tumor","carcinoma", "sarcomas","leukemia", "myeloma","metastasis","seminoma","melanoma")#,"oma")
index <- grep(paste(toMatch,collapse="|"),not_found$sample_attributes, ignore.case = TRUE)

extra_cancer<-not_found[index,c("experiment_acc","sample_acc" )]
saveRDS(extra_cancer, file="/dcs07/hansen/data/recount_ASE/data/cancer_annot_additional.rds")
#------------------------------------------------
#get single cell 
#------------------------------------------------



recount3_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/metadata/Recount3_metadata.tsv")


recount3_metadata2<-as_tibble(recount3_metadata) %>% 
  dplyr::select(study,experiment_acc,external_id,study_abstract,study_title,library_selection,sample_attributes,sample_description)
sra<-split(recount3_metadata2,1:nrow(recount3_metadata2))
names(sra)<-recount3_metadata2$study


potential_single_cell <- c()
for(i in c(1:length(sra))){
  if(grepl("single cell", sra[[i]]$study_abstract[1], ignore.case = T)){
    potential_single_cell <- c(potential_single_cell, sra[[i]]$external_id)
  }
  else{
    if(grepl("RIP-seq", sra[[i]]$study_title[1], ignore.case = T)){
      potential_single_cell <- c(potential_single_cell, sra[[i]]$external_id)
    }
    else{
      if(grepl("ICESeq", sra[[i]]$study_title[1], ignore.case = T)){
        potential_single_cell <- c(potential_single_cell, sra[[i]]$external_id)
      }
      else{
        if(grepl("ICE-seq", sra[[i]]$study_title[1], ignore.case = T)){
          potential_single_cell <- c(potential_single_cell, sra[[i]]$external_id)
        }
    else{
      if(grepl("single-cell", sra[[i]]$study_abstract[1], ignore.case = T)){
        potential_single_cell <- c(potential_single_cell, sra[[i]]$external_id)
      }
      else{
        if(grepl("circRNA", sra[[i]]$study_abstract[1], ignore.case = T)){
          potential_single_cell <- c(potential_single_cell, sra[[i]]$external_id)
        }
        else{
          if(grepl("scRNA", sra[[i]]$study_abstract[1], ignore.case = T)){
            potential_single_cell <- c(potential_single_cell, sra[[i]]$external_id)
          }else{
            if(grepl("single-nucleus", sra[[i]]$study_abstract[1], ignore.case = T)){
              potential_single_cell <- c(potential_single_cell, sra[[i]]$external_id)
            }else{
              if(grepl("single nucleus", sra[[i]]$study_abstract[1], ignore.case = T)){
                potential_single_cell <- c(potential_single_cell, sra[[i]]$external_id)
              }else{
                if(grepl("snRNA", sra[[i]]$study_abstract[1], ignore.case = T)){
                  potential_single_cell <- c(potential_single_cell, sra[[i]]$external_id)
                }else{
                  if(grepl("SEQC", sra[[i]]$study_abstract[1], ignore.case = T)){
                    potential_single_cell <- c(potential_single_cell, sra[[i]]$external_id)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
    }
  }}

saveRDS(potential_single_cell, file = "/dcs07/hansen/data/recount_ASE/data/potential_single_cell.rds")


#---------------------------------------------------------------------------------------------------------------------
#Get tissue parent
#---------------------------------------------------------------------------------------------------------------------

library("rols")

ontology<-readRDS("/dcs07/hansen/data/recount_ASE/data/sra_ontology_term.rds")

uberon <- Ontology("uberon")
uberontrms <- Terms(uberon) ## or Terms("bspo")

df<-termId(uberontrms)
id<-match(ontology$term_id,termId(uberontrms))
id_nona<-id[!is.na(id)]
parent_term<-rep(NA, length(id))
nona_indx<-which(!is.na(id))

for(i in 1:length(id_nona)){
  print(i)
  nn<-id_nona[i]
  kk<-nona_indx[i]
parent_term[kk]<-termLabel(uberontrms[[nn]])
}

saveRDS(parent_term,"~/ASE/tissue.rds")

