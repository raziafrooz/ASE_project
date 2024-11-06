
library(recount3)
setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(ggplot2)


human_projects <- available_projects()
human_projects <- available_projects()
qc_meta<-subset(human_projects, file_source == "gtex" & project_type == "data_sources")
plotFile<-"/dcs07/hansen/data/recount_ASE/metadata/gtex_qc_metadata.csv.gz"


if(!file.exists(plotFile)){


qc_df<-c()

for(ss in qc_meta$project){
  print(ss)
  url<-locate_url(
    ss,
    "data_sources/gtex",
    type = "metadata")
  
  xx1 <- utils::read.delim(file_retrieve(url[1], verbose = FALSE))
  xx3 <- utils::read.delim(file_retrieve(url[3], verbose = FALSE))
  xx4 <- utils::read.delim(file_retrieve(url[4], verbose = FALSE))

  tissue<-full_join(xx1,xx3,by=c("external_id","rail_id"))
  tissue<-full_join(tissue,xx4,by=c("external_id","rail_id"))
  
  qc_df<-rbind(qc_df,tissue)
  rm(tissue)
}

fwrite(qc_df, plotFile)
}else{
  qc_df<-read.csv(plotFile)
}
