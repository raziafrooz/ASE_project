library(data.table)
library(tidyverse)
library(MetBrewer)


ase_df<-fread("/dcs07/hansen/data/recount_ASE/data/sra_ASE_stat.csv.gz")
recount3_metadata<-fread("/dcs07/hansen/data/recount_genotype/new_count_pipeline/new_count_pipeline/metadata/Recount3_metadata.tsv", header= T, sep = "\t",quote="")
# xx<-readRDS("recount_genotype/redo_manuscript_figures/ready_to_plot/all_SRA.annotated.rds")
# saveRDS(xx,"/dcs07/hansen/data/recount_ASE/data/sra_nSNP.rds")
sra_nSNP<-readRDS("/dcs07/hansen/data/recount_ASE/data/sra_nSNP.rds")
sra_nSNP$experiment_acc<-recount3_metadata$experiment_acc[match(sra_nSNP$sample_id,recount3_metadata$external_id)]
ase_df$nSNP<-sra_nSNP$nSNPs[match(ase_df$experiment_acc,sra_nSNP$experiment_acc)]


ase_df$library_layout<-recount3_metadata$library_layout[match(ase_df$experiment_acc,recount3_metadata$experiment_acc)]
ase_df$bc_frag.mode_length<-recount3_metadata$bc_frag.mode_length[match(ase_df$experiment_acc,recount3_metadata$experiment_acc)]
ase_df$star.average_input_read_length<-as.numeric(recount3_metadata$star.average_input_read_length[match(ase_df$experiment_acc,recount3_metadata$experiment_acc)])

ase_df<-ase_df[-which(ase_df$ref_ratio==0),] %>% 
  mutate(overlap= star.average_input_read_length-bc_frag.mode_length)


#Single end should not have bc_frag.mode_length, but some have >0 values--> fix them
ase_df$overlap[which(ase_df$library_layout=="single" & ase_df$bc_frag.mode_length==0)]<-ase_df$star.average_input_read_length[which(ase_df$library_layout=="single" & ase_df$bc_frag.mode_length==0)]

#-----------------------------
#Remove samples that are not bulk and added by mistake and also
#remove cancer samples for now


scRNA<-readRDS("/dcs07/hansen/data/recount_ASE/data/potential_single_cell.rds")
cancer_df<-readRDS("/dcs07/hansen/data/recount_ASE/data/cancer_annot.rds")
cancer_df$external_id<-recount3_metadata$external_id[match(cancer_df[,1],recount3_metadata$sample_acc)]

cancer_df_additional<-readRDS("/dcs07/hansen/data/recount_ASE/data/cancer_annot_additional.rds")
cancer_df_additional$external_id<-recount3_metadata$external_id[match(cancer_df_additional$sample_acc,recount3_metadata$sample_acc)]



bad_samp<-data.frame(external_id=unique(c(cancer_df$external_id,cancer_df_additional$external_id)))
bad_samp_2<-data.frame(external_id=unique(scRNA))


bad_samp$experiment_acc<-recount3_metadata$experiment_acc[match(bad_samp$external_id,recount3_metadata$external_id)]
bad_samp_2$experiment_acc<-recount3_metadata$experiment_acc[match(bad_samp_2$external_id,recount3_metadata$external_id)]

dim(ase_df)
ase_df<-ase_df[-which(ase_df$experiment_acc %in% unique(bad_samp$experiment_acc)),]
dim(ase_df)
ase_df<-ase_df[-which(ase_df$experiment_acc %in% unique(bad_samp_2$experiment_acc)),]


#-----------------------------

single_df <-ase_df %>% filter(library_layout=="single") %>% mutate(het_ratio = nHet/nSNP)

paired_df<-ase_df %>% filter(library_layout!="single") %>% mutate(het_ratio = nHet/nSNP)

nhet_lim_single_df<-800 #remove 10% samples with low nHET=774
nrow(single_df)-sum(single_df$nHet >= nhet_lim_single_df) #4,050 removed
nhet_lim_paired_df<-1000 #quantile(paired_df$nHet, probs = seq(0, 1, by = .025))[[2]] #remove 5% samples with low nHET=1915
nrow(paired_df)-sum(paired_df$nHet >= nhet_lim_paired_df) #2,010 removed



overlap_lim_paired_low<- -100 #quantile(paired_df$overlap, probs = seq(0, 1, by = .05))[[2]] #remove 5% samples with low overlap< -66
overlap_lim_paired_high<- 150 #quantile(paired_df$overlap, probs = seq(0, 1, by = .05))[[20]] #remove 5% samples with high overlap > 130

nrow(paired_df)-sum(paired_df$ref_ratio>= 0.4 &paired_df$ref_ratio<=0.6)

single_filter<-single_df %>%  filter(nHet >= nhet_lim_single_df,
                                     overlap >= 0,
                                     !ref_ratio< 0.4,
                                     !ref_ratio>0.6,
                                     q75 >15)
nrow(single_df)-nrow(single_filter) #5,682 removed --> 5,201 -->5,710

paired_filter<-paired_df %>%  filter(nHet >= nhet_lim_paired_df,
                                     overlap >= overlap_lim_paired_low,
                                     overlap <= overlap_lim_paired_high,
                                     !ref_ratio< 0.4,
                                     !ref_ratio>0.6,
                                     q75 >15)

final_ase<-rbind(paired_filter,single_filter) #remain:106,550 --> 96,163 --> 88,935


nrow(paired_filter) #50,181
nrow(single_filter) #42,800
nrow(final_ase) #92,981



#-----------------------------------------------------------------------------------------------

cell_type<-readRDS("/dcs07/hansen/data/recount_ASE/data/cell_type_annot.rds")
ontology<-readRDS("/dcs07/hansen/data/recount_ASE/data/sra_ontology_term.rds")
parent_term<-readRDS("~/ASE/tissue.rds")

ontology$parent_tissue<-parent_term


final_ase$sample_accession<-recount3_metadata$sample_acc[match(final_ase$experiment_acc,recount3_metadata$experiment_acc)]

#There are multiple DIOD terms for one sample but if the sample is diseased, it should have DOID:4 anyways
disease<-ontology[str_starts(ontology$term_id, "DOID"),]


final_ase$sample_type<-cell_type$sample_type[match(final_ase$sample_accession,cell_type$sample_accession)]
final_ase$disease<-disease$term_id[match(final_ase$sample_accession,disease$sample_accession)]


final_tissue<-final_ase %>% filter(sample_type %in% c("primary cells", "tissue"))
no_annot<-final_ase %>% filter(is.na(sample_type))
no_annot$study<-recount3_metadata$study[match(no_annot$experiment_acc,recount3_metadata$experiment_acc)]

if(!file.exists("/dcs07/hansen/data/recount_ASE/data/toUse_tissue_primaryCell.csv")){
fwrite(final_tissue,"/dcs07/hansen/data/recount_ASE/data/toUse_tissue_primaryCell.csv")
fwrite(no_annot,"/dcs07/hansen/data/recount_ASE/data/toUse_no_annot.csv")

}
