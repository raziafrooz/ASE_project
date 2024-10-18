library(recount3)

sra<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/metadata/all_SRA.csv")
sra[1,]
sample_id<-sra_geno$sample_id[id][k]

human_projects <- available_projects()
i=9
for(i in 1:length(unique(sra$study)) ){
  one_study=sra_geno$study[id][i]
  
proj_info <- subset(
  human_projects,
  project == one_study & project_type == "data_sources"
)

## Create a RSE object at the gene level
rse_gene <- create_rse(proj_info)

## Expand the SRA attributes (see details for more information)
rse_gene <- expand_sra_attributes(rse_gene)

meta_df<-colData(rse_gene)[,c("external_id","study","sra.sample_title","sra.sample_description","sra.sample_attributes",
                     "recount_pred.curated.tissue","recount_pred.pattern.predict.type",
                     "recount_pred.curated.type","recount_pred.pred.type","recount_pred.curated.cell_type",
                     "recount_pred.curated.cell_line","sra_attribute.sample_comment")]


}
colnames(colData(rse_gene))[170:185]
colData(rse_gene)[1,168:175]

projects=data.frame(project=unique(sra_geno$study[id]))
metadata_files <- do.call(rbind, apply(projects, 1, function(x) {
  x <-
    locate_url(
      project = x[["project"]],
      type = "metadata"
    )
  res <- data.frame(t(x))
  colnames(res) <-
    gsub("\\..*", "", gsub("^[a-z]+\\.", "", colnames(res)))
  
  colnames(res)[colnames(res) %in% unique(projects$file_source)] <-
    "project_meta"
  x=res$recount_pred[1]
  res<-read.delim(recount3::file_retrieve(x))
  
  return(res)
}))
save(metadata_files, file="~/test/met.rda")
colnames(metadata_files)
metadata_files[1:2,]
table(metadata_files$curated.cell_type)
