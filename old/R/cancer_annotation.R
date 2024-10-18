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

FirstTenDOBimap <- DOTERM
xx <- as.list(FirstTenDOBimap)


ds_term<-unlist(lapply(xx,function(y) Term(y)))
doid_tb<-tibble(doid=names(xx), temrs=ds_term)

toMatch<-c("cancer", "tumor","carcinoma", "sarcomas","leukemia","oma")
index <- grep(paste(toMatch,collapse="|"),doid_tb$temrs, ignore.case = TRUE)
doid_tb$cancer<-"none_cancer"
doid_tb$cancer[index]<-"cancer"
yy<-unique(doid_tb$doid[which(doid_tb$cancer=="cancer")])


#------------------------------
#Get sample annotation from metaSRA
#------------------------------

#download the sqlite latest version from https://metasra.biostat.wisc.edu/supportpages/download.html
#load it in:
connection<-dbConnect(drv=RSQLite::SQLite(), dbname="/dcs07/hansen/data/recount_ASE/data/metasra.v1-8.sqlite")

rs <- dbSendQuery(connection, "SELECT * FROM mapped_ontology_terms")
rs1 <- dbFetch(rs)
cancer_samp<- rs1 %>% 
  filter(term_id %in% yy) %>%
  dplyr::select(sample_accession) %>% unique()


saveRDS(cancer_samp, file="/dcs07/hansen/data/recount_ASE/data/cancer_annot.rds")

#------------------------------------------------
#get single cell 
#------------------------------------------------



recount3_metadata<-read_tsv("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/metadata/Recount3_metadata.tsv")


recount3_metadata2<-recount3_metadata %>% 
  select(study,experiment_acc,external_id,study_abstract,study_title,library_selection,sample_attributes,sample_description)
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

