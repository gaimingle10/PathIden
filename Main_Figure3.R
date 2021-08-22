### Main Figure3 #####

library(ggplot2)
library(reshape2)
library(ggpubr)
library(patchwork)
library("RColorBrewer")
library(openxlsx)
library(tidyverse)
library('ggbeeswarm')

myfill = brewer.pal(12,"Paired")

dRNA_ratio <- read.xlsx("data/diff_sequencing_ratio.xlsx",sheet = 3)
dRNA_ratio$type <- factor(dRNA_ratio$type,levels = c("BALF","Blood","CSF"))

p1 <- ggplot(dRNA_ratio,aes(type,ratio))+
  geom_beeswarm(aes(color=sequencing),cex = 1.5)+
  scale_color_manual(values = myfill[c(2,6)])+
  labs(x='',y='Ratio of Effective reads')+
  theme_bw()+my_theme+theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60))


dRNA_database <- read.xlsx("data/diff_sequencing_ratio.xlsx",sheet = 2)
dRNA_database.melt <- melt(dRNA_database)
dRNA_database.melt$type <- factor(dRNA_database.melt$type,levels = c("BALF","Blood","CSF"))


p2 <- ggplot(dRNA_database.melt,aes(x=variable,y=log10(value+1e-8))) +
  geom_beeswarm(aes(color=sequencing),cex=1.5)+ #cex用于设置点的密集程度
  scale_color_manual(values = myfill[c(2,6)])+
  labs(x='',y='Log10(Ratio of readsnum)')+
  facet_wrap(.~type,scales = 'free_x')+
  theme_bw()+my_theme+theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60))

plot_grid(p1+theme(legend.position = 'none'), p2, 
          align = 'h', nrow = 1, axis = "l", rel_widths = c(0.5,1.5),
          greedy = TRUE)
