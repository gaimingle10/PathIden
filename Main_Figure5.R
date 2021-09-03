# Started on 2021-4-15
# UPDATE 2021-8-26
# Run as follows:
# Rscript Main_Figure4.R <diff_sequencing_ratio>

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("A minimum of 1 arguments are mandatory: Rscript Main_Figure4.R <diff_sequencing_ratio>", call.=FALSE)
}

### Main Figure5 #####

library(ggplot2)
library(reshape2)
library(ggpubr)
library(patchwork)
library("RColorBrewer")
library(openxlsx)
library(tidyverse)
library('ggbeeswarm')
library(cowplot)

#set theme for ggplot2
my_theme <- theme(axis.ticks.length = unit(0.4,"lines"), 
                  axis.ticks = element_line(color='black'),
                  axis.line = element_line(colour = "black"), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.title = element_text(colour = 'black',size=10,face="bold"),
                  axis.text.y=element_text(colour='black',size=10,face = "bold"),
                  axis.text.x=element_text(colour='black',size=10,face = "bold"),
                  plot.title=element_text(hjust = 0.5))
myfill = brewer.pal(12,"Paired")

diff_sequencing_ratio <- as.character(args[1])
dRNA_ratio <- read.xlsx(diff_sequencing_ratio,sheet = 3)
dRNA_ratio$type <- factor(dRNA_ratio$type,levels = c("BALF","Blood","CSF"))

p1 <- ggplot(dRNA_ratio,aes(type,ratio))+
  geom_beeswarm(aes(color=sequencing),cex = 1.5)+
  scale_color_manual(values = myfill[c(2,6)])+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8))+
  labs(x='',y='Ratio of Effective reads')+
  theme_bw()+my_theme+theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60))


dRNA_database <- read.xlsx("data/diff_sequencing_ratio.xlsx",sheet = 2)
dRNA_database.melt <- melt(dRNA_database)
dRNA_database.melt$type <- factor(dRNA_database.melt$type,levels = c("BALF","Blood","CSF"))


p2 <- ggplot(dRNA_database.melt,aes(x=variable,y=log10(value+1e-8))) +
  geom_beeswarm(aes(color=sequencing),cex=1.5)+ #cex用于设置点的密集程度
  scale_color_manual(values = myfill[c(2,6)],name="Sequence")+
  labs(x='',y='Log10(Ratio of readsnum)')+
  facet_wrap(.~type,scales = 'free_x')+
  theme_bw()+my_theme+theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60))

pdf("Figure5.pdf",height = 6,width = 8)
plot_grid(p1+theme(legend.position = 'none'), p2, 
          align = 'h', nrow = 1, axis = "b", rel_widths = c(0.5,1.5),
          greedy = TRUE,labels = c("A","B"))
dev.off()
print("Done")


