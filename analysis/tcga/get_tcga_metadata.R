library(recount3)
library(data.table)
human_projects <- available_projects()
qc_meta<-subset(human_projects, file_source == "tcga" & project_type == "data_sources")

filePath<-"/dcs07/hansen/data/recount_ASE/data/tcga_recount_metadata.csv"
if(!file.exists(filePath)){
all_metadata<-c()
for(ss in qc_meta$project){
  print(ss)
  url<-locate_url(
    ss,
    "data_sources/tcga",
    type = "metadata")
  xx <-utils::read.delim(file_retrieve(url[1], verbose = FALSE))
  xx2 <-utils::read.delim(file_retrieve(url[3], verbose = FALSE))
  study<-left_join(xx,xx2)
  
  all_metadata<-rbind(all_metadata,study)
}
fwrite(all_metadata,filePath)
}


