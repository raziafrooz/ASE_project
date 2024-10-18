geno_met<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/AggregateFiles/GTEx_Blended_Tissue_Testing_metadata.csv")

geno_met$allGenotypesOutput<-paste0("/dcs04/hansen/data/recount_genotype/pipeline/GTEx_Blended_Tissue_Testing/",geno_met$study, "/GTEX_blended_model/compute_conditional_accuracy/",geno_met$study,"_allGenotypedSamplesAgg.rds" )


write.csv(geno_met, "~/ASE/data/GTEx_geno_metadata.csv")
