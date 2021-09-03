# Started on 2021-4-15
# UPDATE 2021-8-26
# Run as follows:
# Rscript Main_Figure4.R <diff_sequencing_ratio>

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("A minimum of 1 arguments are mandatory: Rscript Main_Figure4.R <diff_sequencing_ratio>", call.=FALSE)
}

### Main Figure4 #####
library(ggplot2)
library(reshape2)
library(ggpubr)
library(patchwork)
library("RColorBrewer")
library(openxlsx)
library(tidyverse)
library(rstatix)

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
diff_ratio <- read.xlsx(diff_sequencing_ratio,sheet = 1)
diff_ratio_cor <- diff_ratio

names(diff_ratio_cor)[c(4,5)] <- c("23s_rRNA","MetaPhlAn")
diff_ratio_cor.melt <- melt(diff_ratio_cor)
names(diff_ratio_cor.melt)[4] <- c("database")

diff_ratio_cor.melt$sequencing[diff_ratio_cor.melt$sequencing=="MetaGenome"] <- "mNGS"
diff_ratio_cor.melt$sequencing[diff_ratio_cor.melt$sequencing=="MetaTranscriptome"] <- "mtNGS"
diff_ratio_cor.melt$sequencing[diff_ratio_cor.melt$sequencing=="ONT cDNA"] <- "mtTGS"
diff_ratio_cor.melt$sequencing <- factor(diff_ratio_cor.melt$sequencing,levels = c("mtTGS","mtNGS","mNGS"))
diff_ratio_cor.melt$sample_type[diff_ratio_cor.melt$sample=="bronchoalveolar lavage fluid sample"] <- "BALF"
diff_ratio_cor.melt$sample_type[diff_ratio_cor.melt$sample=="blood sample"] <- "Blood"
diff_ratio_cor.melt$sample_type[diff_ratio_cor.melt$sample=="cerebrospinal fluid sample"] <- "CSF"

diff_ratio_cor.melt$sample_type <- factor(diff_ratio_cor.melt$sample_type,levels = c("BALF","Blood","CSF"))


diff_ratio_cor.melt$value <- log10(diff_ratio_cor.melt$value+1e-8)
#paired pvalue
stat.test <- diff_ratio_cor.melt %>%
  group_by(database,sample_type) %>%
  wilcox_test(value ~ sequencing,ref.group = "mtTGS") 

stat.test <- stat.test %>%
  add_xy_position(x = "database", dodge = 0.8,step.increase = 0.04)

pdf("Figure4.pdf",height = 6,width = 8)

ggboxplot(
  diff_ratio_cor.melt, x = "database", y = "value", color = "sequencing",
  add = "jitter",facet.by = "sample_type",xlab = "",ylab = "Log10(ratio of reads number)")+
  scale_color_manual(name="Detection method",values = myfill[c(4,6,2)])+
  theme_bw()+my_theme+theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 60))

bxp + 
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0,
    bracket.nudge.y= -0.05,hide.ns = FALSE
  )

dev.off()
print("Done")

