## Main Figure1 ####
library(ggplot2)
library(reshape2)
library(ggpubr)
library(patchwork)
library(reshape2)
library("RColorBrewer")
library(tidyverse)

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

##Figure 1A
ratio <- read.csv('data/ratio_maping.csv',header = T)
names(ratio)[c(2,3)] <- c("DNA","RNA")

ratio_plot <- function(data, plt_title){
  p <- ggplot(data,aes(x=DNA*100,y=RNA*100)) +
    geom_point(stat='identity',color="#1F78B4")+theme_bw()+
    geom_abline(slope=5.6,intercept=0)+
    scale_x_continuous(limits = c(0,12.5))+
    scale_y_continuous(limits = c(0,70))+
    labs(x='Ratio of efficient readsnum in mNGS (%)',
         y='Ratio of efficient readsnum in mtNGS (%)',
         title = plt_title)+my_theme
  return(p)
}
#total
p1 <- ratio_plot(ratio,"Total")

#balf
ratio_alv <- ratio[which(ratio$sample_type=="bronchoalveolar lavage fluid sample"),]
p2 <- ratio_plot(ratio_alv,"BALF")

#blood
ratio_bd <- ratio[which(ratio$sample_type=="blood sample"),]
p3 <- ratio_plot(ratio_bd,"Blood")

#csf
ratio_cfs <- ratio[which(ratio$sample_type=="cerebrospinal fluid sample"),]
p4 <- ratio_plot(ratio_cfs,"CSF")

pA <- p1+p2+p3+p4+plot_layout(nrow=1)
pA

##Figure 1B
clean_readnum <- read.csv('data/clean_readnum.csv',header = T)
meta_data <- read.csv('data/sample_type.csv')
RNA_16s <- read.csv('data/mtNGS/RNA_16s_genus.txt',header = T,sep = '\t',row.names = 1)
DNA_16s <- read.csv('data/mNGS/DNA_16s_genus.txt',header = T,sep = '\t',row.names = 1)

#merge mNGS and mtNGS data with common species
integrate_data1 <- function(RNA_16s,DNA_16s){
  #intersect
  bacteria <- intersect(row.names(DNA_16s),row.names(RNA_16s))
  DNA_16s.in <- DNA_16s[which(row.names(DNA_16s)%in%bacteria),]
  DNA_16s.in$bacteria <- row.names(DNA_16s.in)
  DNA_16s.in.melt <- melt(DNA_16s.in)
  names(DNA_16s.in.melt)[c(2,3)] <- c("sample",'DNA')
  
  RNA_16s.in <- RNA_16s[which(row.names(RNA_16s)%in%bacteria),]
  RNA_16s.in$bacteria <- row.names(RNA_16s.in)
  RNA_16s.in.melt <- melt(RNA_16s.in)
  names(RNA_16s.in.melt)[c(2,3)] <- c("sample",'RNA')
  
  merged_16s <- merge(RNA_16s.in.melt,DNA_16s.in.melt)
  merged_16s$type <- meta_data$sample_type[match(merged_16s$sample,meta_data$name)]
  return(merged_16s)
}
#set colors
cols <- RColorBrewer::brewer.pal(10,'Paired')

#filter data based on average relative abudance
fil_aver <- function(DNA_16s,RNA_16s){
  #calculate relative abundance
  for(i in 1:ncol(DNA_16s)){
    DNA_16s[,i] <- DNA_16s[,i]/clean_readnum$MG[which(clean_readnum$sample==names(DNA_16s)[i])]*100
  }
  for(i in 1:ncol(RNA_16s)){
    RNA_16s[,i] <- RNA_16s[,i]/clean_readnum$MT[which(clean_readnum$sample==names(RNA_16s)[i])]*100
  }
  
  #merge mNGS and mtNGS
  merged_16s <- integrate_data1(RNA_16s,DNA_16s)
  
  #remove zero value both in mNGS and mtNGS
  index <- c()
  for (i in 1:nrow(merged_16s)) {
    if(merged_16s[i,3]!=0 && merged_16s[i,4]!=0) index <- append(index,i)
  }
  merged_16s.f <- merged_16s[index,]
  merged_16s.f$bacteria <- factor(merged_16s.f$bacteria,levels = unique(merged_16s.f$bacteria))
  
  #calculate average Relative abundance and set colors
  bacteria_mean <- merged_16s.f %>%
    group_by(bacteria) %>%
    summarise(mean=mean(DNA+RNA))
  bacteriaTop10 <- bacteria_mean[order(bacteria_mean$mean,decreasing = T)[1:10],]$bacteria
  top10_color <- data.frame(bacteria=bacteriaTop10,color=cols)
  merged_16s.f$cols <- top10_color$color[match(merged_16s.f$bacteria,top10_color$bacteria)]
  merged_16s.f$cols <- as.character(merged_16s.f$cols)
  merged_16s.f$cols[is.na(merged_16s.f$cols)] <- "grey90"
  return(merged_16s.f)
}

#plot function with ggplot2
scatter_plot <- function(data,xlim,ylim,ylabel){
  p <- ggplot(data,aes(x=log10(DNA+1e-8),y=log10(RNA+1e-8))) +
    geom_point(stat='identity',color=data$cols)+theme_bw()+
    geom_abline(slope=1,intercept=0)+
    scale_x_continuous(limits = c(xlim,ylim))+
    scale_y_continuous(limits = c(xlim,ylim))+
    stat_cor(aes(label = paste(..rr.label..,..p.label..,sep = "~`,`~")),label.y=ylabel,color="blue")+
    labs(x='Log10(mNGS)',y='Log10(mtNGS)')+my_theme
  return(p)
}
#total
merged_16s.f <- fil_aver(DNA_16s,RNA_16s)
pB1 <- scatter_plot(merged_16s.f,-6.5,0,0)

#balf
merged_16s.f.balf <- merged_16s.f[which(merged_16s.f$type=="bronchoalveolar lavage fluid sample"),]
pB2 <- scatter_plot(merged_16s.f.balf,-6.5,0,0)

#blood
merged_16s.f.bs <- merged_16s.f[which(merged_16s.f$type=="blood sample"),]
pB3 <- scatter_plot(merged_16s.f.bs,-6.5,0,0)

#cfs
merged_16s.f.cfs <- merged_16s.f[which(merged_16s.f$type=="cerebrospinal fluid sample"),]
pB4 <- scatter_plot(merged_16s.f.cfs,-6.5,0,0)

#construct legends
y=seq(-6.5,-1.5,0.5)
pB5 <- ggplot()+lims(x=c(-1,1))+theme_void()
pB5 <- pB5+annotate("point",x=-1,y=y,shape=16,color=c("grey90",rev(cols)),size=2)+
  annotate("text", x=-0.9,y=y,label=c("Others",rev(as.character(merged_16s.f$bacteria[match(cols,merged_16s.f$cols)]))), hjust=0)

pB <- pB1+pB2+pB3+pB4+pB5+plot_layout(nrow=1)
pB

##Figure 1C
RNA_23s <- read.csv('data/mtNGS/RNA_23s_genus.txt',header = T,sep = '\t',row.names = 1)
DNA_23s <- read.csv('data/mNGS/DNA_23s_genus.txt',header = T,sep = '\t',row.names = 1)

#total
merged_23s.f <- fil_aver(DNA_23s,RNA_23s)
pC1 <- scatter_plot(merged_23s.f,-6.5,0,0)

#balf
merged_23s.f.balf <- merged_23s.f[which(merged_23s.f$type=="bronchoalveolar lavage fluid sample"),]
pC2 <- scatter_plot(merged_23s.f.balf,-6.5,0,0)

#blood
merged_23s.f.bs <- merged_23s.f[which(merged_23s.f$type=="blood sample"),]
pC3 <- scatter_plot(merged_23s.f.bs,-6.5,0,0)

#cfs
merged_23s.f.cfs <- merged_23s.f[which(merged_23s.f$type=="cerebrospinal fluid sample"),]
pC4 <- scatter_plot(merged_23s.f.cfs,-6.5,0,0)

#construct legends
y=seq(-6.5,-1.5,0.5)
pC5 <- ggplot()+lims(x=c(-1,1))+theme_void()
pC5 <- pC5+annotate("point",x=-1,y=y,shape=16,color=c("grey90",rev(cols)),size=2)+
  annotate("text", x=-0.9,y=y,label=c("Others",rev(as.character(merged_23s.f$bacteria[match(cols,merged_23s.f$cols)]))), hjust=0)

pC <- pC1+pC2+pC3+pC4+pC5+plot_layout(nrow=1)
pC

##Figure 1D
RNA_bacteria <- read.table('data/mtNGS/RNA_Bacteria_genus.txt',header = T,row.names = 1)
DNA_bacteria <- read.table('data/mNGS/DNA_Bacteria_genus.txt',header = T,row.names = 1)

#total
merged_bacter.f <- fil_aver(DNA_bacteria,RNA_bacteria)
pD1 <- scatter_plot(merged_bacter.f,-5,1.2,1)

#balf
merged_bacter.f.balf <- merged_bacter.f[which(merged_bacter.f$type=="bronchoalveolar lavage fluid sample"),]
pD2 <- scatter_plot(merged_bacter.f.balf,-5,1.2,1)

#blood
merged_bacter.f.bs <- merged_bacter.f[which(merged_bacter.f$type=="blood sample"),]
pD3 <- scatter_plot(merged_bacter.f.bs,-5,1.2,1)

#cfs
merged_bacter.f.cfs <- merged_bacter.f[which(merged_bacter.f$type=="cerebrospinal fluid sample"),]
pD4 <- scatter_plot(merged_bacter.f.cfs,-5,1.2,1)

#construct legends
y=seq(-6.5,-1.5,0.5)
pD5 <- ggplot()+lims(x=c(-1,1))+theme_void()
pD5 <- pD5+annotate("point",x=-1,y=y,shape=16,color=c("grey90",rev(cols)),size=2)+
  annotate("text", x=-0.9,y=y,label=c("Others",rev(as.character(merged_bacter.f$bacteria[match(cols,merged_bacter.f$cols)]))), hjust=0)

pD <- pD1+pD2+pD3+pD4+pD5+plot_layout(nrow=1)
pD
