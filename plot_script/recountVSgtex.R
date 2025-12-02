
library(data.table)
library(ggplot2)
library(cowplot)
library(tidyverse)
theme_set(theme_cowplot())
theme1= theme(axis.text=element_text(size = 10), axis.title=element_text(size = 20))+
theme_half_open() +
  background_grid()

tissue_abb<-read.table("~/ASE-data/data/gtex_tissue_abbre.txt", sep="\t")

gtex_noWASP<-"/dcs07/hansen/data/recount_ASE/data/gtexVSrecount.csv.gz"
gtex_WASP<-"/dcs07/hansen/data/recount_ASE/data/gtexVSrecount_wasp.csv.gz"


com_plot_no_wasp<-fread(gtex_noWASP)
com_plot_wasp<-fread(gtex_WASP)
colnames(com_plot_wasp)[3:4]<-c("Recount3","GTEx_WASP")
colnames(com_plot_no_wasp)[3:4]<-c("Recount3","GTEx_noWASP")

com_plot<-com_plot_no_wasp %>% full_join(com_plot_wasp) %>% pivot_longer(!c(tissue,sample),values_to="ref_ratio",names_to="pipeline")

com_plot$Tissue<-tissue_abb$V2[match(com_plot$tissue,tissue_abb$V1)]

pdf(file="~/plot/ASE/recountVSgtex2.pdf", width = 10, height = 8)

pp=ggplot(com_plot, aes(y=ref_ratio, x=Tissue))+
  geom_point(alpha=0.2)+
  ylim(c(0.45,0.55))+
  geom_hline(yintercept = 0.5,linetype = "dashed", alpha=0.5, color="red")+theme1+
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(vars(pipeline))+
  labs(y="Median Reference Ratio")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
print(pp)

dev.off()


