library(recount3)
setwd("~/ASE/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(ggplot2)
library(ggridges)

sra_subset<-readRDS("~/test/m.rds")

sra_subset$library_source<-NA
sra_subset$library_selection<-NA
sra_subset$platform_model<-NA
sra_subset$run_published<-NA


for(ss in unique(sra_subset$study)){
  print(ss)
  url<-locate_url(
    ss,
    "data_sources/sra",
    type = "metadata")
  
  xx <-utils::read.delim(file_retrieve(url[1], verbose = FALSE))




  id<-match(sra_subset$sample_id,xx$external_id)
  id<-id[!is.na(id)]
  id2<-match(xx$external_id,sra_subset$sample_id)
  id2<-id2[!is.na(id2)]
  
  
  sra_subset$library_source[id2]<-xx$library_source[id]
  sra_subset$library_selection[id2]<-xx$library_selection[id]
  sra_subset$platform_model[id2]<-xx$platform_model[id]
  sra_subset$run_published[id2]<-xx$run_published[id]
}
sra_subset$run_published<-as_datetime(sra_subset$run_published)
sra_subset$year<-year(sra_subset$run_published)
table(sra_subset$platform_model, sra_subset$year)


#-----------------------------------------------------------------------------------------
#year plot for paired-end sequencing
#-----------------------------------------------------------------------------------------

pdf(file="~/plot/ASE/year_bias_paired.pdf", width = 10, height = 6)

#year vs ref_ratio
p= ggplot(sra_subset[!is.na(sra_subset$library_layout),] %>% filter(library_layout=="paired"),
          aes(x=ref_ratio,y=factor(year), fill=factor(year), color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at paired SRA samples: year vs ref_ratio")
print(p)

#year vs seq length

p= ggplot(sra_subset[!is.na(sra_subset$library_layout),]%>% filter(library_layout=="paired"),
          aes(x=seq_len,y=factor(year), fill=factor(year), color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at paired SRA samples: year vs seq length")
print(p)

#year vs frag len

p= ggplot(sra_subset[!is.na(sra_subset$library_layout),]%>% filter(library_layout=="paired"),
          aes(x=frag_len,y=factor(year), fill=factor(year), color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at paired SRA samples: year vs fragment length")
print(p)


#year vs overlap

p= ggplot(sra_subset[!is.na(sra_subset$library_layout),]%>% filter(library_layout=="paired"),
          aes(x=overlap,y=factor(year), fill=factor(year), color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at paired SRA samples: year vs overlap")
print(p)

p= ggplot(sra_subset[!is.na(sra_subset$library_layout),]%>% filter(library_layout=="paired"),
          aes(x=overlap,y=ref_ratio, color=factor(year)))+
  geom_point(alpha=0.5)+
  labs(title="only looking at paired-end SRA samples: overlap vs ref_ratio each year")
print(p)

dev.off()


#-----------------------------------------------------------------------------------------
#year plot for single-end sequencing
#-----------------------------------------------------------------------------------------



pdf(file="~/plot/ASE/year_bias_single.pdf", width = 10, height = 6)

#year vs ref_ratio

p= ggplot(sra_subset[!is.na(sra_subset$library_layout),] %>% filter(library_layout=="single"),
          aes(x=ref_ratio,y=factor(year), fill=factor(year), color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at single-end SRA samples: year vs ref_ratio")
print(p)

#year vs seq length


p= ggplot(sra_subset[!is.na(sra_subset$library_layout),]%>% filter(library_layout=="single"),
          aes(x=seq_len,y=factor(year), fill=factor(year), color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at single-end SRA samples: year vs seq length")
print(p)


p= ggplot(sra_subset[!is.na(sra_subset$library_layout),]%>% filter(library_layout=="single"),
          aes(x=seq_len,y=ref_ratio, color=factor(year)))+
  geom_point(alpha=0.5)+
  labs(title="only looking at single-end SRA samples: year vs seq length")
print(p)


dev.off()



#-----------------------------------------------------------------------------------------
#platform plot for paired-end sequencing
#-----------------------------------------------------------------------------------------


pdf(file="~/plot/ASE/platform_bias_paired.pdf", width = 10, height = 6)


#platform vs ref_ratio

p= ggplot(sra_subset[!is.na(sra_subset$library_layout),] %>% filter(library_layout=="paired"),
          aes(x=ref_ratio,y=platform_model, fill=platform_model, color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at paired-end SRA samples: platform_model vs ref_ratio")
print(p)

#platform_model vs seq length

p= ggplot(sra_subset[!is.na(sra_subset$library_layout),] %>% filter(library_layout=="paired"),
          aes(x=seq_len,y=platform_model, fill=platform_model, color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at paired-end SRA samples: platform_model vs sequencing length")
print(p)
#platform_model vs frag len
p= ggplot(sra_subset[!is.na(sra_subset$library_layout),] %>% filter(library_layout=="paired"),
          aes(x=frag_len,y=platform_model, fill=platform_model, color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at paired-end SRA samples: platform_model vs fragment len")
print(p)

#platform_model vs overlap
p= ggplot(sra_subset[!is.na(sra_subset$library_layout),] %>% filter(library_layout=="paired"),
          aes(x=overlap,y=platform_model, fill=platform_model, color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at paired-end SRA samples: platform_model vs overlap")
print(p)

#platform_model overlap vs ref_ratio
p= ggplot(sra_subset[!is.na(sra_subset$library_layout),]%>% filter(library_layout=="paired"),
          aes(x=overlap,y=ref_ratio, color=platform_model))+
  geom_point(alpha=0.5)+
  labs(title="only looking at paired-end SRA samples: overlap vs ref_ratio each platform_model")
print(p)



# p= ggplot(sra_subset[!is.na(sra_subset$library_layout),]%>% filter(library_layout=="paired") %>% mutate(aFC=abs(0.5-ref_ratio)),
#           aes(y=platform_model,x=aFC,color=platform_model))+
#   geom_boxplot(alpha=0.5)+
#   geom_jitter(aes(color=platform_model),alpha=0.4, width=0.2, size=2)
#   labs(title="only looking at paired-end SRA samples: platform_model vs aFC ")
# print(p)

dev.off()


#-----------------------------------------------------------------------------------------
#platform plot for single-end sequencing
#-----------------------------------------------------------------------------------------


pdf(file="~/plot/ASE/platform_bias_single.pdf", width = 10, height = 6)

#platform vs ref_ratio

p= ggplot(sra_subset[!is.na(sra_subset$library_layout),] %>% filter(library_layout=="single"),
          aes(x=ref_ratio,y=platform_model, fill=platform_model, color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at single-end SRA samples: platform_model vs ref_ratio")
print(p)

#platform_model vs seq length

p= ggplot(sra_subset[!is.na(sra_subset$library_layout),] %>% filter(library_layout=="single"),
          aes(x=seq_len,y=platform_model, fill=platform_model, color=library_layout))+
  geom_density_ridges(alpha=0.5,jittered_points = TRUE)+
  labs(title="only looking at single-end SRA samples: platform_model vs sequencing length")
print(p)

#platform_model seq_len vs ref_ratio
p= ggplot(sra_subset[!is.na(sra_subset$library_layout),]%>% filter(library_layout=="single"),
          aes(x=seq_len,y=ref_ratio, color=platform_model))+
  geom_point(alpha=0.5)+
  labs(title="only looking at paired-end SRA samples: seq_len vs ref_ratio each platform_model")
print(p)


dev.off()