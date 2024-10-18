setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(MetBrewer)


tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("~/ASE-data/data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

tissues_names$file_name<-paste0("/dcl01/hansen/data/arazi/ASE/dbGap/GTEx_Analysis_v8_ASE_counts_by_tissue/", tissues_names$file_name)
wasp<-tissues_names


tissues_names<-as.data.frame(list.files(path ="/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/"))
colnames(tissues_names)<-"file_name"
tissues_names$abb<-sapply(strsplit(gsub("\\.","-",tissues_names$file_name),"-"), function(xx){ xx[2]  })

tissue_abb<-read.table("~/ASE-data/data/gtex_tissue_abbre.txt", sep="\t")
colnames(tissue_abb)<-c("name","abb")
tissues_names$full_name<-tissue_abb$name[match(tissues_names$abb,tissue_abb$abb)]

tissues_names$file_name<-paste0("/dcl01/hansen/data/gtex_ase/GTEx_Analysis_v8_ASE_counts_by_tissue/", tissues_names$file_name)
no_wasp<-tissues_names
rm(tissues_names)

plotFile<-"~/test/waspVSno_wasp.csv.gz"
com_plot<-c()
tissues<-wasp$full_name
for(ss in 2:length(tissues)){
  study<-tissues[ss]
  print(study)
  
  wasp_1<- fread(wasp$file_name[wasp$full_name==study][1])
  colnames(wasp_1)[1:2]<- c("chr", "start")
  
  
  no_wasp_1<- fread(no_wasp$file_name[no_wasp$full_name==study][1])
  colnames(no_wasp_1)[1:2]<- c("chr", "start")
  
  
 
    
   
    join_ase<- inner_join(wasp_1,no_wasp_1, by=c("chr", "start","SAMPLE_ID"))
    rm(no_wasp_1)
    rm(wasp_1)
    
    com_plot2<- join_ase %>%
      group_by(SAMPLE_ID) %>% 
      summarize(wasp=median(REF_RATIO.x,na.rm=T),
             no_wasp=median(REF_RATIO.y,na.rm=T))
      
    
    com_plot<-rbind(com_plot,com_plot2)
    rm(com_plot2)
  }
fwrite(com_plot,plotFile)





plotFile_no_wasp<-fread("/dcs07/hansen/data/recount_ASE/data/gtexVSrecount.csv.gz")
plotFile_wasp<-fread("/dcs07/hansen/data/recount_ASE/data/gtexVSrecount_wasp.csv.gz")
colnames(plotFile_no_wasp)[3:4]<-c("no_wap_recount","no_wasp_gtex")
colnames(plotFile_wasp)[3:4]<-c("wap_recount","wasp_gtex")

plot_df<-inner_join(plotFile_wasp,plotFile_no_wasp)
plot_df<-plot_df %>%
  pivot_longer(!c(tissue,sample), names_to = "pipeline", values_to = "ref_ratio")




pdf(file="~/plot/ASE/waspVSnowasp.pdf", width = 10, height = 6)


ggplot(plot_df, aes(x=pipeline, y=ref_ratio))+
  geom_boxplot()+
  labs(title="Boxplot of median ref ratio of samples in GTEx with and without WASP correction",
       subtitle="recount was subseted to common SNPs with GTEx in each condition")
dev.off()








