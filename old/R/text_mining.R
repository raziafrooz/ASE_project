

setwd("~/ASE/")
library(tidyverse)
library(ggplot2)
library(ggridges)
library(data.table)


recount3_metadata<-read_tsv("/dcs04/hansen/data/recount_genotype/PCA/SRA/Recount3_metadata.tsv")


recount3_metadata2<-recount3_metadata %>% select(study,external_id,study_abstract,study_title,library_selection,sample_attributes,sample_description)
sra<-split(recount3_metadata2,1:nrow(recount3_metadata2))
names(sra)<-recount3_metadata$study
rm(recount3_metadata)
library_selection <- lapply(sra, function(isra){
  isra$library_selection
})

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
sra[["SRP063998"]] <- NULL
sra <- sra[!names(sra) %in% potential_single_cell]
#saveRDS(potential_single_cell, file = "data/potential_single_cell.rds")
# sra <- lapply(sra, function(isra){
#   keep <- isra$library_selection == "cDNA" | isra$sra.library_selection == "RT-PCR" | isra$sra.library_selection == "PCR"
#   isra[ ,keep]
# })

# nsra <- lapply(sra, function(rse){
#   dim(rse)[2]
# })
# 
# nsra <- unlist(nsra)
# sra <- sra[nsra >= 20]

sra_metadata <- lapply(sra, function(rse){
  rse$sample_attributes
}) 
sra_metadata[["SRP014574"]] <- sra[["SRP014574"]]$sample_description
mtD_df_list <- list()
for(n in c(1:length(sra_metadata))){
  print(n)
  mtD <- sra_metadata[[n]]
  
  if(sum(is.na(mtD)) == length(mtD)){
    mtD_df_list[[n]] <- NA
  } else{
    mtD <- strsplit(mtD, split = "\\|")
    lmtD <- lapply(mtD, function(imt){ length(imt)})
    lmtD <- unlist(lmtD)
    colnames_df <- c()
    not_zero <- which.max(lmtD)
    mtD_df <- data.frame(matrix(NA, nrow = length(mtD), ncol = length(mtD[[not_zero]])))
    for(i in c(1:length(mtD[[not_zero]]))){
      colnames_df[i] <- unlist(strsplit(mtD[[not_zero]][i], split = ";;"))[1]
    }
    colnames(mtD_df) <- colnames_df
    mtD <- lapply(mtD, function(imtD){
      imtD <- strsplit(imtD, ";;")
      a <- lapply(imtD, function(iimtD){
        data_vec <- iimtD[2]
        names(data_vec) <- iimtD[1]
        data_vec
      })
      unlist(a)
    })
    for(j in c(1:length(mtD))){
      if(length(mtD[[j]]) > 0){
        for(k in colnames(mtD_df)){
          mtD_df[j, k] <- mtD[[j]][k]
        }}else{
          mtD_df[j, ] <- NA
        }
    }
    mtD_df_list[[n]] <- mtD_df}
}

names(mtD_df_list) <- names(sra)
metaData_cols <- lapply(mtD_df_list, function(imt){
  colnames(imt)
})

metaData_cols <- unlist(metaData_cols)
metaData_cols <- unique(metaData_cols)
metaData_cols <- data.frame(metaData_cols)

diseases_if_found <- list()
for(i in c(1:length(mtD_df_list))){
  imt <- mtD_df_list[[i]]
  index <- grep("disease", colnames(imt), ignore.case = TRUE)
  if(length(index) > 0){
    diseases_if_found[[i]] <- imt[ ,index]
  }else{
    diseases_if_found[[i]] <- NA
  }
}
names(diseases_if_found) <- names(mtD_df_list)
diseases_if_found[["ERP010930"]] <- diseases_if_found[["ERP010930"]]$disease
diseases_if_found[["ERP013700"]] <- diseases_if_found[["ERP013700"]]$disease
diseases_if_found[["ERP104602"]] <- diseases_if_found[["ERP104602"]]$disease
diseases_if_found[["ERP106610"]] <- diseases_if_found[["ERP106610"]]$disease
diseases_if_found[["ERP109255"]] <- diseases_if_found[["ERP109255"]]$disease
diseases_if_found[["ERP111913"]] <- diseases_if_found[["ERP111913"]]$disease
diseases_if_found[["SRP052896"]] <- diseases_if_found[["SRP052896"]]$disease_stage
diseases_if_found[["SRP069976"]] <- diseases_if_found[["SRP069976"]]$disease
diseases_if_found[["SRP055438"]] <- diseases_if_found[["SRP055438"]]$disease
diseases_if_found[["SRP075592"]] <- diseases_if_found[["SRP075592"]]$disease
diseases_if_found[["SRP078515"]] <- diseases_if_found[["SRP078515"]]$disease
diseases_if_found[["SRP091981"]] <- diseases_if_found[["SRP091981"]]$disease
diseases_if_found[["SRP148497"]] <- diseases_if_found[["SRP148497"]]$disease
diseases_if_found[["SRP174668"]] <- diseases_if_found[["SRP174668"]]$disease
diseases_if_found[["SRP113755"]] <- diseases_if_found[["SRP113755"]]$disease
diseases_if_found[["SRP117312"]] <- diseases_if_found[["SRP117312"]]$disease

for(i in c(1:length(diseases_if_found))){
  if(is.na(diseases_if_found[[i]])){
    index <- grep("cancer", colnames(mtD_df_list[[i]]), ignore.case = TRUE)
    if(length(index) > 0){
      diseases_if_found[[i]] <- mtD_df_list[[i]][ ,index]}
  }
}

for(i in c(1:length(diseases_if_found))){
  if(is.na(diseases_if_found[[i]])){
    index <- grep("diagnosis", colnames(mtD_df_list[[i]]), ignore.case = TRUE)
    if(length(index) > 0){
      diseases_if_found[[i]] <- mtD_df_list[[i]][ ,index]}
  }
}
diseases_if_found[["SRP042228"]] <- diseases_if_found[["SRP042228"]]$diagnosis
diseases_if_found[["SRP048801"]] <- diseases_if_found[["SRP048801"]]$diagnosis
diseases_if_found[["SRP072494"]] <- diseases_if_found[["SRP072494"]]$diagnosis
diseases_if_found[["SRP113470"]] <- diseases_if_found[["SRP113470"]]$diagnosis
diseases_if_found[["SRP096757"]] <- diseases_if_found[["SRP096757"]]$diagnosis
diseases_if_found[["SRP129004"]] <- diseases_if_found[["SRP129004"]]$diagnosis

for(i in c(1:length(diseases_if_found))){
  if(is.na(diseases_if_found[[i]])){
    index <- grep("histology", colnames(mtD_df_list[[i]]), ignore.case = TRUE)
    if(length(index) > 0){
      diseases_if_found[[i]] <- mtD_df_list[[i]][ ,index]}
  }
}

for(i in c(1:length(diseases_if_found))){
  if(is.na(diseases_if_found[[i]])){
    index <- grep("pathology", colnames(mtD_df_list[[i]]), ignore.case = TRUE)
    if(length(index) > 0){
      diseases_if_found[[i]] <- mtD_df_list[[i]][ ,index]}
  }
}

for(i in c(1:length(diseases_if_found))){
  if(is.na(diseases_if_found[[i]])){
    index <- grep("tumor", colnames(mtD_df_list[[i]]), ignore.case = TRUE)
    if(length(index) > 0){
      diseases_if_found[[i]] <- mtD_df_list[[i]][ ,index]}
  }
}

diseases_if_found[["SRP157215"]] <- diseases_if_found[["SRP157215"]]$`tumor type`


source_if_found <- lapply(mtD_df_list, function(imt){
  index <- grep("source", colnames(imt), ignore.case = TRUE)
  if(length(index) > 0){
    imt[ ,index]
  }else{
    NA
  }
})

for(i in c(1:length(source_if_found))){
  if(is.na(source_if_found[[i]])){
    index <- grep("origin", colnames(mtD_df_list[[i]]), ignore.case = TRUE)
    if(length(index) > 0){
      source_if_found[[i]] <- mtD_df_list[[i]][ ,index]}
  }
}

tissue_if_found <- lapply(mtD_df_list, function(imt){
  index <- grep("tissue", colnames(imt), ignore.case = TRUE)
  if(length(index) > 0){
    imt[ ,index]
  }else{
    NA
  }
})

for(i in c(1:length(tissue_if_found))){
  if(is.na(tissue_if_found[[i]])){
    index <- grep("body", colnames(mtD_df_list[[i]]), ignore.case = TRUE)
    index <- index[!index %in% c("body mass index", "body type")]
    if(length(index) > 0){
      tissue_if_found[[i]] <- mtD_df_list[[i]][ ,index]}
  }
}

cell_if_found <- lapply(mtD_df_list, function(imt){
  index <- grep("cell", colnames(imt), ignore.case = TRUE)
  if(length(index) > 0){
    imt[ ,index]
  }else{
    NA
  }
})

biopsy_if_found <- lapply(mtD_df_list, function(imt){
  index <- grep("biopsy", colnames(imt), ignore.case = TRUE)
  if(length(index) > 0){
    imt[ ,index]
  }else{
    NA
  }
})

organ_if_found <- lapply(mtD_df_list, function(imt){
  index <- grep("organ", colnames(imt), ignore.case = TRUE)
  index <- index[!index %in% c("organism", "organism status", "organisms")]
  if(length(index) > 0){
    imt[ ,index]
  }else{
    NA
  }
})

description_if_found <- lapply(mtD_df_list, function(imt){
  index <- grep("Description", colnames(imt), ignore.case = TRUE)
  if(length(index) > 0){
    imt[ ,index]
  }else{
    NA
  }
})

# Manually fix the studies with multiple columns
tissue_if_found[["DRP000987"]] <- rep("peripheral blood", dim(sra[["DRP000987"]])[2])
tissue_if_found[["DRP001953"]] <- rep("peripheral blood", dim(sra[["DRP001953"]])[2])
tissue_if_found[["ERP000546"]] <- c(mtD_df_list[["ERP000546"]]$OrganismPart[1:16], mtD_df_list[["ERP000546"]]$`SRA accession`[17:dim(mtD_df_list[["ERP000546"]])[1]])
tissue_if_found[["ERP006077"]] <- mtD_df_list[["ERP006077"]]$Description
tissue_if_found[["ERP007185"]] <- mtD_df_list[["ERP007185"]]$Description
tissue_if_found[["ERP010003"]] <- rep("blood", dim(mtD_df_list[["ERP010003"]])[2])
tissue_if_found[["SRP215939"]] <- rep("colon", dim(mtD_df_list[["SRP215939"]])[1])
tissue_if_found[["SRP198979"]] <- mtD_df_list[["SRP198979"]]$condition
tissue_if_found[["SRP183071"]] <- rep("Diffuse large B cell lymphoma", dim(mtD_df_list[["SRP183071"]])[1])
tissue_if_found[["SRP182842"]] <- mtD_df_list[["SRP182842"]]$tisse
tissue_if_found[["ERP004043"]] <- rep("kidney cancer", dim(mtD_df_list[["ERP004043"]])[1])
tissue_if_found[["ERP004043"]][colData(sra[["ERP004043"]])$sra.experiment_acc == "ERX324194"] <- "kidney"
tissue_if_found[["ERP010795"]] <- mtD_df_list[["ERP010795"]]$Alias
tissue_if_found[["ERP010795"]][grep("EGC", tissue_if_found[["ERP010795"]])] <- "gastric cancer"
tissue_if_found[["ERP010795"]][grep("HGD", tissue_if_found[["ERP010795"]])] <- "gastric dyplasia (high)"
tissue_if_found[["ERP010795"]][grep("LGD", tissue_if_found[["ERP010795"]])] <- "gastric dyplasia (low)"
tissue_if_found[["ERP010930"]] <- mtD_df_list[["ERP010930"]]$`ENA checklist`
tissue_if_found[["ERP010930"]][!mtD_df_list[["ERP010930"]]$disease == "Grade 2"] <- mtD_df_list[["ERP010930"]]$disease[!mtD_df_list[["ERP010930"]]$disease == "Grade 2"]
tissue_if_found[["ERP010930"]][37] <- "oligoastrocytoma"
tissue_if_found[["ERP012979"]] <- rep("thyroid", dim(mtD_df_list[["ERP012979"]])[1])
tissue_if_found[["SRP212704"]] <- rep("prostrate", dim(mtD_df_list[["SRP212704"]])[1])
tissue_if_found[["SRP212704"]][grep("T", sra[["SRP212704"]]$sra.experiment_title)] <- "prostrate cancer"
tissue_if_found[["SRP191427"]] <- rep("caudate nucleus", dim(sra[["SRP191427"]])[2])
tissue_if_found[["SRP189540"]] <- mtD_df_list[["SRP189540"]]$isolate
tissue_if_found[["SRP185703"]] <- rep("gastrointestinal stromal tumors", dim(sra[["SRP185703"]])[2])
tissue_if_found[["SRP174991"]][tissue_if_found[["SRP174991"]] == "non-tumor"] <- "liver"
tissue_if_found[["SRP174991"]][tissue_if_found[["SRP174991"]] == "HCC"] <- "liver cancer"
tissue_if_found[["SRP171066"]] <- rep("SSC-40 cells", length(tissue_if_found[["SRP171066"]]))
tissue_if_found[["SRP164655"]] <- rep("breast cancer", length(tissue_if_found[["SRP164655"]]))
tissue_if_found[["SRP154573"]][is.na(tissue_if_found[["SRP154573"]])] <- "glioma tumor"
tissue_if_found[["SRP154388"]] <- mtD_df_list[["SRP154388"]]$`body site`
tissue_if_found[["SRP149433"]] <- rep("breast cancer primary-derived xenograft", length(tissue_if_found[["SRP149433"]]))
tissue_if_found[["ERP013498"]] <- rep("glioblastoma", dim(sra[["ERP013498"]])[2])
tissue_if_found[["ERP014531"]] <- rep("monocytes", dim(sra[["ERP014531"]])[2])
tissue_if_found[["ERP015139"]] <- mtD_df_list[["ERP015139"]]$Description
tissue_if_found[["ERP020523"]] <- mtD_df_list[["ERP020523"]]$`ArrayExpress-CellType`
tissue_if_found[["ERP020535"]] <- rep("duodenum", dim(sra[["ERP020535"]])[2])
tissue_if_found[["ERP022034"]] <- mtD_df_list[["ERP022034"]]$disease
tissue_if_found[["ERP022368"]] <- rep("human monocyte derived macrophages", dim(sra[["ERP022368"]])[2])
tissue_if_found[["ERP023321"]] <- rep("prostrate cancer", dim(sra[["ERP023321"]])[2])
tissue_if_found[["ERP105366"]] <- mtD_df_list[["ERP105366"]]$Title
tissue_if_found[["ERP105366"]][grep("ALL", tissue_if_found[["ERP105366"]])] <- "acute lymphocytic leukemia"
tissue_if_found[["SRP175707"]] <- rep("prostrate cancer xenografts", dim(sra[["SRP175707"]])[2])
tissue_if_found[["SRP131768"]] <- rep("human foreskin fibroblasts", dim(sra[["SRP131768"]])[2])
tissue_if_found[["SRP014830"]] <- rep("breast cancer", dim(sra[["SRP014830"]])[2])
tissue_if_found[["SRP035524"]] <- mtD_df_list[["SRP035524"]]$body_site
tissue_if_found[["SRP042186"]][tissue_if_found[["SRP042186"]] == "WAT"] <- "white adipose tissue"
tissue_if_found[["SRP042186"]][tissue_if_found[["SRP042186"]] == "MPC"] <- "mesenchymal progenitor cells"
tissue_if_found[["SRP042186"]][tissue_if_found[["SRP042186"]] == "BAT"] <- "brown adipose tissue"
tissue_if_found[["SRP119923"]] <- mtD_df_list[["SRP119923"]]$isolate
tissue_if_found[["SRP076716"]] <- mtD_df_list[["SRP076716"]]$`cell or tumor type`
tissue_if_found[["SRP131119"]] <- sra[["SRP131119"]]$sra.experiment_title
tissue_if_found[["SRP131119"]][grep("hMSC", tissue_if_found[["SRP131119"]])] <- "hMSC"
tissue_if_found[["SRP131119"]][grep("A673", tissue_if_found[["SRP131119"]])] <- "A673"
tissue_if_found[["SRP131119"]][grep("ES2", tissue_if_found[["SRP131119"]])] <- "ES2"
tissue_if_found[["SRP131119"]][grep("SK-ES", tissue_if_found[["SRP131119"]])] <- "SK-ES"
tissue_if_found[["SRP131083"]] <- sra[["SRP131083"]]$sra.experiment_title
tissue_if_found[["SRP131083"]][grep("blood", tissue_if_found[["SRP131083"]])] <- "blood"
tissue_if_found[["SRP131083"]][grep("nasal", tissue_if_found[["SRP131083"]])] <- "nasal epithelium"
tissue_if_found[["SRP099844"]][mtD_df_list[["SRP099844"]]$tissue == "Polyp"] <- "Duodenum polyp"
tissue_if_found[["SRP099053"]][tissue_if_found[["SRP099053"]] == "Tumor tissue"] <- "hepatocellular carcinoma"
tissue_if_found[["SRP099053"]][tissue_if_found[["SRP099053"]] == "Non-neoplastic liver tissue"] <- "liver"
tissue_if_found[["SRP051472"]] <- sra[["SRP051472"]]$sra.experiment_title
tissue_if_found[["SRP051472"]][grep("undiff", tissue_if_found[["SRP051472"]])] <- "hESC"
tissue_if_found[["SRP051472"]][grep("meso", tissue_if_found[["SRP051472"]])] <- "mesoderm"
tissue_if_found[["SRP051472"]][grep("cardio", tissue_if_found[["SRP051472"]])] <- "cardiomyocyte"
tissue_if_found[["SRP051472"]][grep("progen", tissue_if_found[["SRP051472"]])] <- "cardiac progenitor"
tissue_if_found[["SRP051249"]] <- sra[["SRP051249"]]$sra.experiment_title
tissue_if_found[["SRP051249"]][grep("Stomach", tissue_if_found[["SRP051249"]])] <- "stomach"
tissue_if_found[["SRP051249"]][grep("Heart", tissue_if_found[["SRP051249"]])] <- "heart"
tissue_if_found[["SRP051249"]][grep("Lung", tissue_if_found[["SRP051249"]])] <- "lung"
tissue_if_found[["SRP051249"]][grep("Kidney", tissue_if_found[["SRP051249"]])] <- "kidney"
tissue_if_found[["SRP051249"]][grep("Intestine", tissue_if_found[["SRP051249"]])] <- "intestine"
tissue_if_found[["SRP051249"]][grep("Adrenal", tissue_if_found[["SRP051249"]])] <- "adrenal"
tissue_if_found[["SRP028530"]] <- sra[["SRP028530"]]$sra.experiment_title
tissue_if_found[["SRP028530"]][grep("PEO1 cell line", tissue_if_found[["SRP028530"]])] <- "PEO1 cell line"
tissue_if_found[["SRP028530"]][grep("C4-4 cell line", tissue_if_found[["SRP028530"]])] <- "ovarian cancer"
tissue_if_found[["SRP018008"]] <- sra[["SRP018008"]]$sra.experiment_title
tissue_if_found[["SRP018008"]][grep("Hepatocellular Carcinoma", tissue_if_found[["SRP018008"]])] <- "hepatocellular carcinoma"
tissue_if_found[["SRP018008"]][grep("bladder cancer", tissue_if_found[["SRP018008"]])] <- "bladder cancer"
tissue_if_found[["SRP017262"]] <- rep("acute myeloid leukemia", dim(sra[["SRP017262"]])[2])
tissue_if_found[["SRP014574"]] <- rep("gastric cancer cell line", dim(sra[["SRP014574"]])[2])
tissue_if_found[["SRP011085"]] <- sra[["SRP011085"]]$sra.sample_title
tissue_if_found[["SRP011085"]][grep("_h", tissue_if_found[["SRP011085"]])] <- "HUVEC"
tissue_if_found[["SRP011085"]][grep("_o", tissue_if_found[["SRP011085"]])] <- "OKF6"
tissue_if_found[["SRP011085"]][grep("co", tissue_if_found[["SRP011085"]])] <- "OKF6"
tissue_if_found[["SRP011085"]][grep("ch", tissue_if_found[["SRP011085"]])] <- "HUVEC"
tissue_if_found[["ERP108822"]] <- rep("breast cancer PDX", dim(sra[["ERP108822"]])[2])
tissue_if_found[["ERP107715"]] <- mtD_df_list[["ERP107715"]]$`organism part`
tissue_if_found[["ERP106876"]] <- mtD_df_list[["ERP106876"]]$Title
tissue_if_found[["ERP106876"]][grep("MOLM-16", tissue_if_found[["ERP106876"]])] <- "acute myeloid leukemia"
tissue_if_found[["ERP106876"]][grep("LP-1", tissue_if_found[["ERP106876"]])] <- "multiple myeloma"
tissue_if_found[["ERP106876"]][grep("EOL-1", tissue_if_found[["ERP106876"]])] <- "chronic eosinophilic leukemia"
tissue_if_found[["ERP105597"]] <- mtD_df_list[["ERP105597"]]$Alias
CNS_tumor <- c("Sample02_rna", "Sample03_rna", "Sample05_rna", "Sample08_rna", "Sample13_rna", "Sample15_rna", "Sample16_rna", "Sample23_rna", "Sample24_rna_1","Sample34_rna", "Sample35_rna_2", "Sample38_rna", "Sample41_rna", "Sample46_rna", "Sample50_rna", "Sample51_rna", "Sample22_rna", "Sample24_rna_2", "Sample28_rna")
tissue_if_found[["ERP105597"]][tissue_if_found[["ERP105597"]] %in% CNS_tumor] <- "CNS tumor"
tissue_if_found[["ERP105597"]][tissue_if_found[["ERP105597"]] == "Sample20_rna"] <- "lymphoma"
tissue_if_found[["ERP105597"]][grep("Sample", tissue_if_found[["ERP105597"]])] <- "extracranial tumor"
tissue_if_found[["ERP015593"]] <- rep("EVSA-T, LnCAP, HCC-70 and ZR-75-1 cell lines", dim(sra[["ERP015593"]])[2])
tissue_if_found[["ERP010372"]] <- rep("human bronchial epithelial cells", dim(sra[["ERP010372"]])[2])
tissue_if_found[["ERP010372"]][grep("CF", mtD_df_list[["ERP010372"]]$Title)] <- paste(tissue_if_found[["ERP010372"]][grep("CF", mtD_df_list[["ERP010372"]]$Title)], "cystic fibrosis")
tissue_if_found[["SRP133789"]] <- paste(mtD_df_list[["SRP133789"]]$histology, mtD_df_list[["SRP133789"]]$source_name)
tissue_if_found[["SRP158639"]][tissue_if_found[["SRP158639"]]== "adjacent tissue of cancer"] <- "pancreas"
tissue_if_found[["SRP069212"]][tissue_if_found[["SRP069212"]]== "adjacent normal"] <- "liver"
tissue_if_found[["SRP075759"]] <- paste(mtD_df_list[["SRP075759"]]$`hbl status`, mtD_df_list[["SRP075759"]]$source_name)
tissue_if_found[["SRP113619"]] <- mtD_df_list[["SRP113619"]]$`brain region`
tissue_if_found[["SRP073789"]] <- paste(mtD_df_list[["SRP073789"]]$`biopsy site`, mtD_df_list[["SRP073789"]]$`normal/tumor`)
tissue_if_found[["SRP058243"]] <- paste(mtD_df_list[["SRP058243"]]$source_name, mtD_df_list[["SRP058243"]]$tissue)
tissue_if_found[["SRP133891"]][tissue_if_found[["SRP133891"]]== "Turmor"] <- "gastric cancer"
tissue_if_found[["SRP133891"]][tissue_if_found[["SRP133891"]]== "Normal"] <- "stomach"
tissue_if_found[["SRP045225"]][tissue_if_found[["SRP045225"]]== "healthy"] <- "lung"
tissue_if_found[["SRP113755"]] <- paste(mtD_df_list[["SRP113755"]]$disease, mtD_df_list[["SRP113755"]]$tissue)
tissue_if_found[["SRP171137"]] <- paste("Schneiderian papillomas", tissue_if_found[["SRP171137"]])
tissue_if_found[["SRP186168"]] <- mtD_df_list[["SRP186168"]]$isolate
tissue_if_found[["SRP068976"]] <- paste(mtD_df_list[["SRP068976"]]$diagnosis, mtD_df_list[["SRP068976"]]$tissue)
tissue_if_found[["SRP056696"]] <- paste(mtD_df_list[["SRP056696"]]$`tumor status`, mtD_df_list[["SRP056696"]]$tissue)
diseases_if_found[["SRP059172"]] <- mtD_df_list[["SRP059172"]]$`history of brucellosis`
cell_if_found[["SRP059172"]] <- NA
tissue_if_found[["ERP109255"]] <- mtD_df_list[["ERP109255"]]$`organism part`

# Fixing labels
diseases_if_found[["SRP191427"]] <- sra[["SRP191427"]]$sra.experiment_title
diseases_if_found[["SRP216012"]] <- sra[["SRP216012"]]$sra.experiment_title
diseases_if_found[["SRP136057"]] <- sra[["SRP136057"]]$sra.experiment_title
diseases_if_found[["SRP145493"]] <- paste("PTSD", mtD_df_list[["SRP145493"]]$PTSD)
diseases_if_found[["SRP155483"]] <-mtD_df_list[["SRP155483"]]$phenotype
diseases_if_found[["SRP034732"]] <- mtD_df_list[["SRP034732"]]$treatment
diseases_if_found[["SRP040547"]]  <-mtD_df_list[["SRP040547"]]$treatment
diseases_if_found[["SRP052896"]] <- mtD_df_list[["SRP052896"]]$disease
diseases_if_found[["SRP053296"]] <- mtD_df_list[["SRP053296"]]$`mi type`
diseases_if_found[["SRP056733"]] <- mtD_df_list[["SRP056733"]]$agent
diseases_if_found[["SRP056840"]] <- mtD_df_list[["SRP056840"]]$`sirs outcomes`
diseases_if_found[["SRP191427"]] <- sra[["SRP191427"]]$sra.experiment_title
diseases_if_found[["SRP216012"]] <- sra[["SRP216012"]]$sra.experiment_title
diseases_if_found[["SRP136057"]] <- sra[["SRP136057"]]$sra.experiment_title
diseases_if_found[["SRP145493"]] <- paste("PTSD", mtD_df_list[["SRP145493"]]$PTSD)
diseases_if_found[["SRP155483"]] <- mtD_df_list[["SRP155483"]]$phenotype
diseases_if_found[["SRP034732"]] <- mtD_df_list[["SRP034732"]]$treatment
diseases_if_found[["SRP040547"]]  <- mtD_df_list[["SRP040547"]]$treatment
diseases_if_found[["SRP049340"]]  <- mtD_df_list[["SRP049340"]]$treatment
diseases_if_found[["SRP051848"]]  <- mtD_df_list[["SRP051848"]]$condition
diseases_if_found[["SRP052896"]] <- mtD_df_list[["SRP052896"]]$disease
diseases_if_found[["SRP053101"]] <- paste(mtD_df_list[["SRP053101"]]$time, mtD_df_list[["SRP053101"]]$`type of surgical procedure`)
diseases_if_found[["SRP056477"]] <- mtD_df_list[["SRP056477"]]$genotype
diseases_if_found[["SRP056822"]] <- mtD_df_list[["SRP056822"]]$group
diseases_if_found[["SRP059039"]] <- paste(mtD_df_list[["SRP059039"]]$group, mtD_df_list[["SRP059039"]]$organisms)
diseases_if_found[["SRP063006"]] <- mtD_df_list[["SRP063006"]]$`hpv status`
source_if_found[["SRP063006"]] <- mtD_df_list[["SRP063006"]]$`tumour site`
diseases_if_found[["SRP063875"]] <- mtD_df_list[["SRP063875"]]$activated
diseases_if_found[["SRP089857"]] <- mtD_df_list[["SRP089857"]]$subject
diseases_if_found[["SRP117629"]] <- mtD_df_list[["SRP117629"]]$treatment
diseases_if_found[["SRP120040"]] <- mtD_df_list[["SRP120040"]]$subject
diseases_if_found[["SRP092544"]] <- paste("Ebola", mtD_df_list[["SRP092544"]]$isolate)
diseases_if_found[["SRP082973"]] <- mtD_df_list[["SRP082973"]]$`smoking status`
diseases_if_found[["SRP077046"]] <- mtD_df_list[["SRP077046"]]$`clinical condition`
diseases_if_found[["SRP077975"]]  <- mtD_df_list[["SRP077975"]]$`clinical information`
diseases_if_found[["SRP078152"]] <- mtD_df_list[["SRP078152"]]$treatment
diseases_if_found[["SRP074736"]] <- mtD_df_list[["SRP074736"]]$condition
diseases_if_found[["SRP072769"]] <- mtD_df_list[["SRP072769"]]$condition
diseases_if_found[["SRP071758"]] <- paste(mtD_df_list[["SRP071758"]]$`copd status`, mtD_df_list[["SRP071758"]]$`dysplasia status`)
diseases_if_found[["SRP064353"]] <- mtD_df_list[["SRP064353"]]$`overexpression of`
diseases_if_found[["SRP064515"]] <- mtD_df_list[["SRP064515"]]$infection
diseases_if_found[["SRP058181"]] <- sra[["SRP058181"]]$sra.experiment_title
diseases_if_found[["SRP062278"]] <- mtD_df_list[["SRP062278"]]$treatment
diseases_if_found[["SRP065865"]] <- mtD_df_list[["SRP065865"]]$condition
diseases_if_found[["SRP078450"]] <- mtD_df_list[["SRP078450"]]$`infection status`
diseases_if_found[["SRP078912"]] <- mtD_df_list[["SRP078912"]]$condition
diseases_if_found[["SRP083918"]] <- paste(mtD_df_list[["SRP083918"]]$group, mtD_df_list[["SRP083918"]]$outcome, mtD_df_list[["SRP083918"]]$timepoint)
diseases_if_found[["SRP091886"]] <- paste(mtD_df_list[["SRP091886"]]$infection, mtD_df_list[["SRP091886"]]$`time point`)
diseases_if_found[["SRP094553"]] <- mtD_df_list[["SRP094553"]]$`functional status/outcome`
diseases_if_found[["SRP094100"]] <- paste(mtD_df_list[["SRP094100"]]$shRNA, mtD_df_list[["SRP094100"]]$treatment)
diseases_if_found[["SRP098715"]] <- sra[["SRP098715"]]$sra.experiment_title
diseases_if_found[["SRP098758"]] <- mtD_df_list[["SRP098758"]]$group
diseases_if_found[["SRP099111"]] <- sra[["SRP099111"]]$sra.experiment_title
diseases_if_found[["SRP101726"]] <- mtD_df_list[["SRP101726"]]$condition
diseases_if_found[["SRP101737"]] <- mtD_df_list[["SRP101737"]]$condition
diseases_if_found[["SRP158730"]] <- mtD_df_list[["SRP158730"]]$`histologic subtype`
diseases_if_found[["SRP136108"]] <- paste(mtD_df_list[["SRP136108"]]$exposure, mtD_df_list[["SRP136108"]]$time_in_hours, "h")
diseases_if_found[["SRP126560"]] <- mtD_df_list[["SRP126560"]]$`genotype/variation`
diseases_if_found[["SRP114864"]] <- sra[["SRP114864"]]$sra.experiment_title
diseases_if_found[["SRP113245"]] <- mtD_df_list[["SRP113245"]]$condition
diseases_if_found[["SRP100787"]] <- mtD_df_list[["SRP100787"]]$condition
diseases_if_found[["SRP100947"]] <- mtD_df_list[["SRP100947"]]$`cultured with`
diseases_if_found[["SRP188219"]] <- sra[["SRP188219"]]$sra.experiment_title
diseases_if_found[["SRP188219"]][grep("AF_RA", diseases_if_found[["SRP188219"]])] <- "Atrial fibrillation right atrium"
diseases_if_found[["SRP188219"]][grep("AF_LA", diseases_if_found[["SRP188219"]])] <- "Atrial fibrillation left atrium"
diseases_if_found[["SRP188219"]][grep("SR_LA", diseases_if_found[["SRP188219"]])] <- "Sinus rhythm left atrium"
diseases_if_found[["SRP188219"]][grep("SR_RA", diseases_if_found[["SRP188219"]])] <- "Sinus rhythm right atrium"
diseases_if_found[["SRP189239"]] <- mtD_df_list[["SRP189239"]]$isolate
diseases_if_found[["SRP151000"]] <- mtD_df_list[["SRP151000"]]$`yu2 strain infection`
diseases_if_found[["SRP136694"]] <- mtD_df_list[["SRP136694"]]$infection
diseases_if_found[["SRP201347"]] <- mtD_df_list[["SRP201347"]]$isolate
diseases_if_found[["SRP116272"]] <- mtD_df_list[["SRP116272"]]$group
diseases_if_found[["SRP162654"]] <- mtD_df_list[["SRP162654"]]$treatment
diseases_if_found[["SRP197353"]] <- mtD_df_list[["SRP197353"]]$`nafld activity score`
diseases_if_found[["SRP197353"]] <- mtD_df_list[["SRP197353"]]$`nafld activity score`
diseases_if_found[["SRP197353"]][!diseases_if_found[["SRP197353"]] == "0"] <- "non alcoholic fatty liver disease"
diseases_if_found[["SRP197353"]][diseases_if_found[["SRP197353"]] == "0"] <- "normal"
diseases_if_found[["SRP149925"]] <- mtD_df_list[["SRP149925"]]$agent
diseases_if_found[["SRP155574"]] <- paste("TGFb replacement", mtD_df_list[["SRP155574"]]$`replacement of tgfb`)
diseases_if_found[["SRP170013"]] <- mtD_df_list[["SRP170013"]]$treatment
diseases_if_found[["SRP149723"]] <- mtD_df_list[["SRP149723"]]$genotype
diseases_if_found[["SRP179696"]] <- mtD_df_list[["SRP179696"]]$isolate
diseases_if_found[["SRP191103"]] <- paste(mtD_df_list[["SRP191103"]]$bmi, mtD_df_list[["SRP191103"]]$`gastric motility`)
diseases_if_found[["SRP119561"]] <- mtD_df_list[["SRP119561"]]$`study group`
diseases_if_found[["SRP162411"]] <- mtD_df_list[["SRP162411"]]$treatment
diseases_if_found[["SRP155182"]] <- sra[["SRP155182"]]$sra.experiment_title
diseases_if_found[["SRP140557"]] <- paste(mtD_df_list[["SRP140557"]]$group, mtD_df_list[["SRP140557"]]$visit)
diseases_if_found[["SRP086612"]] <- paste(mtD_df_list[["SRP086612"]]$tetracycline, mtD_df_list[["SRP086612"]]$treatment)
diseases_if_found[["SRP086626"]] <- mtD_df_list[["SRP086626"]]$treatment
diseases_if_found[["SRP074274"]] <- mtD_df_list[["SRP074274"]]$infection
diseases_if_found[["SRP076235"]] <- sra[["SRP076235"]]$sra.experiment_title
organ_if_found[["SRP078262"]] <- mtD_df_list[["SRP078262"]]$`fetal region`
diseases_if_found[["SRP092402"]] <- paste(mtD_df_list[["SRP092402"]]$`disease state`, mtD_df_list[["SRP092402"]]$treatmentresult)
diseases_if_found[["SRP114315"]] <- mtD_df_list[["SRP114315"]]$isolate
diseases_if_found[["SRP118922"]] <- mtD_df_list[["SRP118922"]]$genotype
diseases_if_found[["SRP120376"]] <- sra[["SRP120376"]]$sra.experiment_title
diseases_if_found[["SRP123014"]] <- mtD_df_list[["SRP123014"]]$isolate
diseases_if_found[["SRP123014"]][grep("S", diseases_if_found[["SRP123014"]])] <- "normal"
diseases_if_found[["SRP123014"]][grep("T", diseases_if_found[["SRP123014"]])] <- "heroin dependent"
diseases_if_found[["SRP123611"]] <- sra[["SRP123611"]]$sra.experiment_title
diseases_if_found[["SRP123611"]][grep("Control", diseases_if_found[["SRP123611"]])] <- "normal"
diseases_if_found[["SRP123611"]][grep("Aneurysm", diseases_if_found[["SRP123611"]])] <- "aneurysm"
diseases_if_found[["SRP131083"]] <- sra[["SRP131083"]]$sra.experiment_title
diseases_if_found[["SRP131083"]][grep("case", diseases_if_found[["SRP131083"]])] <- "asthma"
diseases_if_found[["SRP131083"]][grep("control", diseases_if_found[["SRP131083"]])] <- "normal"
diseases_if_found[["SRP136102"]] <- mtD_df_list[["SRP136102"]]$phenotype
diseases_if_found[["SRP220377"]] <- mtD_df_list[["SRP220377"]]$`disease state`
diseases_if_found[["SRP212956"]] <- mtD_df_list[["SRP212956"]]$`patient group`
diseases_if_found[["SRP217937"]] <- mtD_df_list[["SRP217937"]]$group
diseases_if_found[["SRP192714"]] <- paste(mtD_df_list[["SRP192714"]]$exposure, mtD_df_list[["SRP192714"]]$time)
diseases_if_found[["SRP156762"]] <- mtD_df_list[["SRP156762"]]$`viral infection`
diseases_if_found[["SRP158994"]] <- paste("HCV", mtD_df_list[["SRP158994"]]$`hcv viremia`, mtD_df_list[["SRP158994"]]$hcvgroup, mtD_df_list[["SRP158994"]]$`time point`)
diseases_if_found[["SRP140558"]] <- paste(mtD_df_list[["SRP140558"]]$group, mtD_df_list[["SRP140558"]]$visit)
diseases_if_found[["SRP144647"]] <- mtD_df_list[["SRP144647"]]$allergy_status
diseases_if_found[["SRP159247"]] <- mtD_df_list[["SRP159247"]]$treatment
diseases_if_found[["SRP160881"]] <- mtD_df_list[["SRP160881"]]$sampletype
diseases_if_found[["SRP160881"]][grep("Patient", diseases_if_found[["SRP160881"]])] <- "STAT2 mutation"
diseases_if_found[["SRP162023"]] <- paste("HantavaxTM", mtD_df_list[["SRP162023"]]$`vaccination administration`)
diseases_if_found[["SRP162414"]] <- paste(mtD_df_list[["SRP162414"]]$treatment, mtD_df_list[["SRP162414"]]$hours_post_treatment)
diseases_if_found[["SRP162492"]] <- paste("acute rejection at 3m", mtD_df_list[["SRP162492"]]$`acr at 3m`)
diseases_if_found[["SRP162411"]] <- paste(mtD_df_list[["SRP162411"]]$hours_post_treatment, mtD_df_list[["SRP162411"]]$treatment)
diseases_if_found[["SRP162500"]] <- paste(mtD_df_list[["SRP162500"]]$hours_post_treatment, mtD_df_list[["SRP162500"]]$treatment)
diseases_if_found[["SRP165679"]] <- sra[["SRP165679"]]$sra.experiment_title
diseases_if_found[["SRP179648"]] <- paste(mtD_df_list[["SRP179648"]]$condition, mtD_df_list[["SRP179648"]]$illumination)
diseases_if_found[["SRP188301"]]<- paste("live attenuated influenza vaccine", mtD_df_list[["SRP188301"]]$`vaccination timepoint`)
diseases_if_found[["SRP200575"]]<- paste(mtD_df_list[["SRP200575"]]$treatment, mtD_df_list[["SRP200575"]]$`time point`)
diseases_if_found[["SRP220377"]] <- mtD_df_list[["SRP220377"]]$`disease state`
diseases_if_found[["SRP198242"]] <- paste(sra[["SRP198242"]]$sra.experiment_title, mtD_df_list[["SRP198242"]]$adjuvant)
diseases_if_found[["SRP220377"]] <- mtD_df_list[["SRP220377"]]$`disease state`
diseases_if_found[["SRP151577"]] <- mtD_df_list[["SRP151577"]]$treatment
diseases_if_found[["SRP135788"]] <- mtD_df_list[["SRP135788"]]$condition
diseases_if_found[["SRP098758"]] <- mtD_df_list[["SRP098758"]]$group
diseases_if_found[["SRP040472"]] <- sra[["SRP040472"]]$sra.sample_name 
diseases_if_found[["SRP040070"]] <- sra[["SRP040070"]]$sra.sample_name 
diseases_if_found[["SRP009316"]] <- sra[["SRP009316"]]$sra.experiment_title
diseases_if_found[["SRP009316"]][grep("cell line", diseases_if_found[["SRP009316"]])] <- "long term Burkitt lymphoma"
diseases_if_found[["SRP009316"]][grep("SLN", diseases_if_found[["SRP009316"]])] <- "sporadic Burkitt Lymphoma"
diseases_if_found[["ERP107762"]] <- sra[["ERP107762"]]$sra.sample_title
diseases_if_found[["ERP016347"]] <- mtD_df_list[["ERP016347"]]$`developmental stage`
cell_if_found[["ERP016347"]] <- sra[["ERP016347"]]$sra.sample_attributes
cell_if_found[["ERP016347"]][grep("chondrocyte", cell_if_found[["ERP016347"]])] <- "chondrocytes"
cell_if_found[["ERP016347"]][grep("osteocyte", cell_if_found[["ERP016347"]])] <- "osteocytes"
cell_if_found[["ERP016347"]][grep("mesenchymal stem cell", cell_if_found[["ERP016347"]])] <- "mesenchymal stem cells"
cell_if_found[["ERP016347"]][grep("tenocyte", cell_if_found[["ERP016347"]])] <- "tenocytes"
cell_if_found[["ERP012527"]] <- sra[["ERP012527"]]$sra.sample_acc.x
cell_if_found[["ERP012527"]][grep("ERS874706", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874708", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874710", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874712", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874714", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874716", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874718", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874720", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874722", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874724", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874726", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874728", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874730", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874732", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874734", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874736", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874738", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874740", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874742", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874744", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874746", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874748", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874750", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS874752", cell_if_found[["ERP012527"]])] <- "BPH1"
cell_if_found[["ERP012527"]][grep("ERS", cell_if_found[["ERP012527"]])] <- "RWPE1"
diseases_if_found[["ERP011275"]] <- sra[["ERP011275"]]$sra.sample_title
diseases_if_found[["ERP011275"]][grep("mock", diseases_if_found[["ERP011275"]])] <- "normal"
diseases_if_found[["ERP011275"]][grep("WA314dYopM_", diseases_if_found[["ERP011275"]])] <- "WA314 YopM deleted"
diseases_if_found[["ERP011275"]][grep("WA314_", diseases_if_found[["ERP011275"]])] <- "WA314"
diseases_if_found[["SRP096589"]] <- sra[["SRP096589"]]$sra.experiment_title
diseases_if_found[["SRP078156"]] <- sra[["SRP078156"]]$sra.experiment_title
diseases_if_found[["SRP098758"]] <- mtD_df_list[["SRP098758"]]$group
tissue_if_found[["ERP008967"]] <- sra[["ERP008967"]]$sra.sample_description
diseases_if_found[["SRP098758"]] <- mtD_df_list[["SRP098758"]]$group
tissue_if_found[["ERP010227"]] <- sra[["ERP010227"]]$sra.sample_title
tissue_if_found[["ERP010227"]][grep("501", tissue_if_found[["ERP010227"]])] <- "airway epithelium"
tissue_if_found[["ERP010227"]][grep("401", tissue_if_found[["ERP010227"]])] <- "airway epithelium"
tissue_if_found[["ERP010227"]][grep("1001", tissue_if_found[["ERP010227"]])] <- "airway epithelium"
tissue_if_found[["ERP010227"]][grep("E208", tissue_if_found[["ERP010227"]])] <- "hepatocytes"
tissue_if_found[["ERP010227"]][grep("E313", tissue_if_found[["ERP010227"]])] <- "hepatocytes"
tissue_if_found[["ERP010227"]][grep("E269", tissue_if_found[["ERP010227"]])] <- "hepatocytes"
tissue_if_found[["ERP010227"]][grep("E321", tissue_if_found[["ERP010227"]])] <- "hepatocytes"
diseases_if_found[["ERP010227"]] <- sra[["ERP010227"]]$sra.sample_title
diseases_if_found[["ERP006215"]] <- sra[["ERP006215"]]$sra.sample_title

i <- names(sra)[1]
study <- sra[[i]]$study
run <- sra[[i]]$external_id
# sample <- colnames(sra[[i]])
tissue <- tissue_if_found[[i]]
biopsy <- biopsy_if_found[[i]]
cell <- cell_if_found[[i]]
organ <- organ_if_found[[i]]
disease <- diseases_if_found[[i]]
source <- source_if_found[[i]]
des <- description_if_found[[i]]
tissue_df <- data.frame(cbind(study,run, tissue, organ, biopsy, cell, disease, source, des))
for(i in names(sra)[197460:length(sra)]){
  print(i)
  study <- sra[[i]]$study
  run <- sra[[i]]$external_id
  # sample <- colnames(sra[[i]])
  tissue <- tissue_if_found[[i]]
  biopsy <- biopsy_if_found[[i]]
  cell <- cell_if_found[[i]]
  organ <- organ_if_found[[i]]
  disease <- diseases_if_found[[i]]
  source <- source_if_found[[i]]
  des <- description_if_found[[i]]
  if(is.data.frame(tissue)){
    tissue <- do.call(paste, tissue)
  }
  if(is.data.frame(des)){
    des <- do.call(paste, des)
  }
  if(is.data.frame(cell)){
    cell <- do.call(paste, cell)
  }
  if(is.data.frame(organ)){
    organ <- do.call(paste, organ)
  }
  if(is.data.frame(disease)){
    disease <- do.call(paste, disease)
  }
  if(is.data.frame(source)){
    source <- do.call(paste, source)
  }
  if(is.data.frame(biopsy)){
    biopsy <- do.call(paste, biopsy)
  }
  tmp <- data.frame(cbind(study,run, tissue, organ, biopsy, cell, disease, source, des))
  tissue_df <- rbind(tissue_df, tmp)
}

write.csv(tissue_df, "~/plot/ASE/recount3_metadata1.csv")


dim(tissue_df)
