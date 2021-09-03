# Started on 2021-4-15
# UPDATE 2021-8-26
# Run as follows:
# Rscript Main_Figure1-3.R <ratio_maping> <clean_readnum> <sample_type> <RNA_16s_genus> <DNA_16s_genus>
# <RNA_23s_genus> <DNA_23s_genus> <RNA_Bacteria_genus> <DNA_Bacteria_genus> <RNA_virus_species> <DNA_virus_species>
# <RNA_18s_species> <DNA_18s_species> <RNA_ARGs_reads> <DNA_ARGs_reads> <aro_index>

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 16) {
  stop("A minimum of 16 arguments are mandatory: Rscript Main_Figure1-3.R <ratio_maping> <clean_readnum> 
       <sample_type> <RNA_16s_genus> <DNA_16s_genus> <RNA_23s_genus> <DNA_23s_genus> 
       <RNA_Bacteria_genus> <DNA_Bacteria_genus> <RNA_virus_species> <DNA_virus_species> 
       <RNA_18s_species> <DNA_18s_species> <RNA_ARGs_reads> <DNA_ARGs_reads> ", call.=FALSE)
}

## Main Figure1 ####
library(ggplot2)
library(reshape2)
library(ggpubr)
library(patchwork)
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
ratio_mapping <- as.character(args[1])
ratio <- read.csv(ratio_mapping,header = T)
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

pdf('Figure1_A.pdf',height = 3.5, width = 16)
p1+p2+p3+p4+plot_layout(nrow=1)
dev.off()
print("Completed figure1 A")

##Figure 1B
clean_readnum <- as.character(args[2])
sample_type <- as.character(args[3])
RNA_16s_genus <- as.character(args[4])
DNA_16s_genus <- as.character(args[5])

clean_readnum <- read.csv(clean_readnum,header = T)
meta_data <- read.csv(sample_type)
RNA_16s <- read.csv(RNA_16s_genus,header = T,sep = '\t',row.names = 1)
DNA_16s <- read.csv(DNA_16s_genus,header = T,sep = '\t',row.names = 1)

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

pdf('Figure1_B.pdf',height = 3.5, width = 16)
pB1+pB2+pB3+pB4+pB5+plot_layout(nrow=1)
dev.off()
print("Completed figure1 B")

##Figure 1C
RNA_23s_genus <- as.character(args[6])
DNA_23s_genus <- as.character(args[7])

RNA_23s <- read.csv(RNA_23s_genus,header = T,sep = '\t',row.names = 1)
DNA_23s <- read.csv(DNA_23s_genus,header = T,sep = '\t',row.names = 1)

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

pdf('Figure1_C.pdf',height = 3.5, width = 16)
pC1+pC2+pC3+pC4+pC5+plot_layout(nrow=1)
dev.off()
print("Completed figure1 C")

##Figure 1D
RNA_Bacteria_genus <- as.character(args[8])
DNA_Bacteria_genus <- as.character(args[9])

RNA_bacteria <- read.table(RNA_Bacteria_genus,header = T,row.names = 1)
DNA_bacteria <- read.table(DNA_Bacteria_genus,header = T,row.names = 1)

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

pdf('Figure1_D.pdf',height = 3.5, width = 16)
pD1+pD2+pD3+pD4+pD5+plot_layout(nrow=1)
dev.off()
print("Completed figure1 D")

#############Main Figure2#######
#Figure 2A
RNA_virus_species <- as.character(args[10])
DNA_virus_species <- as.character(args[11])

RNA_virus <- read.csv(RNA_virus_species,header = T,sep = '\t',row.names = 1)
DNA_virus <- read.csv(DNA_virus_species,header = T,sep = '\t',row.names = 1)

#total
merged_virus.f <- fil_aver(DNA_virus,RNA_virus)
pA1 <- scatter_plot(merged_virus.f,-6.5,1,1)

#balf
merged_virus.f.balf <- merged_virus.f[which(merged_virus.f$type=="bronchoalveolar lavage fluid sample"),]
pA2 <- scatter_plot(merged_virus.f.balf,-6.5,1,1)

#blood
merged_virus.f.bs <- merged_virus.f[which(merged_virus.f$type=="blood sample"),]
pA3 <- scatter_plot(merged_virus.f.bs,-6.5,1,1)

#cfs
merged_virus.f.cfs <- merged_virus.f[which(merged_virus.f$type=="cerebrospinal fluid sample"),]
pA4 <- scatter_plot(merged_virus.f.cfs,-6.5,1,1)

#construct legends
y=seq(-6.5,0,0.6)
pA5 <- ggplot()+lims(x=c(-1,1))+theme_void()
pA5 <- pA5+annotate("point",x=-1,y=y,shape=16,color=c("grey90",rev(cols)),size=2)+
  annotate("text", x=-0.9,y=y,label=c("Others",rev(as.character(merged_virus.f$bacteria[match(cols,merged_virus.f$cols)]))), hjust=0)

pdf('Figure2_A.pdf',height = 3.5, width = 16)
pA1+pA2+pA3+pA4+pA5+plot_layout(nrow=1)
dev.off()
print("Completed figure2 A")


# plot Venn
extract_data <- function(DNA_16s){
  DNA_16s.t <- data.frame(DNA_16s)
  DNA_16s.t$spcies <- row.names(DNA_16s.t)
  DNA_16s.tm <- melt(DNA_16s.t)
  return(DNA_16s.tm)
}
DNA_16s_data <- extract_data(DNA_16s)
RNA_16s_data <- extract_data(RNA_16s)
DNA_23s_data <- extract_data(DNA_23s)
RNA_23s_data <- extract_data(RNA_23s)
DNA_bacteria_data <- extract_data(DNA_bacteria)
RNA_bacteria_data <- extract_data(RNA_bacteria)
DNA_virus_data <- extract_data(DNA_virus)
RNA_virus_data <- extract_data(RNA_virus)

library(VennDiagram)
myfill = brewer.pal(12,"Paired")
venn_plot <- function(RNA,DNA){
  venn.plot <- venn.diagram(
    x = list(RNA,DNA),
    category.names = c("mtNGS","mNGS"),
    filename = NULL,
    
    #output features
    imagetype = "png",
    width = 480,
    height = 480,
    resolution = 300,
    compression = "lzw",
    
    #Circles
    lwd = 2,
    lty = "blank",
    fill = myfill[c(5,1)],
    cex = .6, #size
    fontfamily = "sans",
    fontface = "bold", 
    
    #set names
    cat.cex = .6, #size
    cat.fontface = "bold", 
    cat.default.pos = "outer",
    cat.pos = c(-27,27)
  )
  return(venn.plot)
}

#extract species name
extract_name <- function(DNA_virus_data,title){
  DNA_virus_data.f <- DNA_virus_data[which(DNA_virus_data$value!=0),]
  sample_type <- meta_data$name[meta_data$sample_type==title]
  name <- unique(DNA_virus_data.f$spcies[which(DNA_virus_data.f$variable %in% sample_type)])
  return(name)
}

#covert data type
library(cowplot)
library(ggplotify)
glist2ggplot <- function(venn.plot){
  vb = as_grob(venn.plot)
  gb = as.ggplot(vb)
  return(gb)
}


#venn between mNGS vs. mtNGS
venn_compare <- function(RNA_16s,DNA_16s,RNA_16s_data,DNA_16s_data){
  #total
  total.venn <- venn_plot(row.names(RNA_16s),row.names(DNA_16s))
  #BALF
  DNA_name <- extract_name(DNA_16s_data,"bronchoalveolar lavage fluid sample")
  RNA_name <- extract_name(RNA_16s_data,"bronchoalveolar lavage fluid sample")
  BALF.venn <- venn_plot(RNA_name,DNA_name)
  #CSF
  DNA_name <- extract_name(DNA_16s_data,"cerebrospinal fluid sample")
  RNA_name <- extract_name(RNA_16s_data,"cerebrospinal fluid sample")
  CFS.venn <- venn_plot(RNA_name,DNA_name)
  #Blood
  DNA_name <- extract_name(DNA_16s_data,"blood sample")
  RNA_name <- extract_name(RNA_16s_data,"blood sample")
  Blood.venn <- venn_plot(RNA_name,DNA_name)
  
  total.venn.plot <- glist2ggplot(total.venn)
  BALF.venn.plot <- glist2ggplot(BALF.venn)
  CFS.venn.plot <- glist2ggplot(CFS.venn)
  Blood.venn.plot <- glist2ggplot(Blood.venn)
  venn16s <- total.venn.plot+BALF.venn.plot+Blood.venn.plot+CFS.venn.plot+plot_layout(nrow = 1)
  return(venn16s)
}

##FigureS3
venn16s <- venn_compare(RNA_16s,DNA_16s,RNA_16s_data,DNA_16s_data)
venn23s <- venn_compare(RNA_23s,DNA_23s,RNA_23s_data,DNA_23s_data)
vennbacteria <- venn_compare(RNA_bacteria,DNA_bacteria,RNA_bacteria_data,DNA_bacteria_data)

pdf("Supplementary_figure3.pdf")
venn16s/venn23s/vennbacteria
dev.off()


##Figure 2B
pdf('Figure2_B.pdf',height = 3.5, width = 12)
venn_compare(RNA_virus,DNA_virus,RNA_virus_data,DNA_virus_data)
dev.off()
print("Completed figure2 B")

##Figure 2C
x <- intersect(merged_virus.f.balf$bacteria,merged_virus.f.bs$bacteria)
y <- intersect(x,merged_virus.f.cfs$bacteria)
virus_data <- merged_virus.f[merged_virus.f$bacteria %in% y,]

#top20
virus_mean <- virus_data %>%
  group_by(bacteria) %>%
  summarise(mean=mean(DNA+RNA))
virusTop20 <- virus_mean[order(virus_mean$mean,decreasing = T)[1:20],]$bacteria
virus_data <- virus_data[which(virus_data$bacteria %in% virusTop20),]
virus_data.melt <- melt(virus_data)
virus_data.melt$sequencing[virus_data.melt$variable=="DNA"] <- "mNGS"
virus_data.melt$sequencing[virus_data.melt$variable=="RNA"] <- "mtNGS"
virus_data.melt$sequencing[virus_data.melt$variable=="cDNA"] <- "ONT cDNA"
virus_data.melt$sequencing <- factor(virus_data.melt$sequencing,levels = c("mNGS","mtNGS","ONT cDNA"))

virus_data.melt$sample_type[virus_data.melt$type=="bronchoalveolar lavage fluid sample"] <- "BALF"
virus_data.melt$sample_type[virus_data.melt$type=="blood sample"] <- "Blood"
virus_data.melt$sample_type[virus_data.melt$type=="cerebrospinal fluid sample"] <- "CFS"
virus_data.melt$sample_type <- factor(virus_data.melt$sample_type,
                                      levels = c("BALF","Blood","CFS"))

pdf('Figure2_C.pdf',height = 8, width = 10)
ggplot(virus_data.melt,aes(bacteria,log10(value+1e-8),color=sequencing))+
  geom_boxplot(position = "dodge2")+geom_jitter(position = position_jitterdodge(0.2))+
  scale_color_manual(values = myfill[c(1,5)])+
  facet_wrap(.~sample_type,scales = 'free_x')+
  coord_flip()+labs(x='',y="Log10(Relative Abundance)")+
  theme_bw()+my_theme
dev.off()
print("Completed figure2 C")

#########Main Figure3#############
##Figrue 3A
RNA_18s_species <- as.character(args[12])
DNA_18s_species <- as.character(args[13])

RNA_18s <- read.csv(RNA_18s_species,header = T,sep = '\t',row.names = 1)
DNA_18s <- read.csv(DNA_18s_species,header = T,sep = '\t',row.names = 1)

DNA_18s_data <- extract_data(DNA_18s)
RNA_18s_data <- extract_data(RNA_18s)

pdf("Figure3_A.pdf",height = 3.5,width = 12)
venn_compare(RNA_18s,DNA_18s,RNA_18s_data,DNA_18s_data)
dev.off()
print("Completed figure3 A")

##Figure 3B
merged_18s.f <- fil_aver(DNA_18s,RNA_18s)
pB1 <- scatter_plot(merged_18s.f,-6.5,1,1)

#balf
merged_18s.f.balf <- merged_18s.f[which(merged_18s.f$type=="bronchoalveolar lavage fluid sample"),]
pB2 <- scatter_plot(merged_18s.f.balf,-6.5,1,1)

#blood
merged_18s.f.bs <- merged_18s.f[which(merged_18s.f$type=="blood sample"),]
pB3 <- scatter_plot(merged_18s.f.bs,-6.5,1,1)

#cfs
merged_18s.f.cfs <- merged_18s.f[which(merged_18s.f$type=="cerebrospinal fluid sample"),]
pB4 <- scatter_plot(merged_18s.f.cfs,-6.5,1,1)

#construct legends
y=seq(-6.5,0,0.6)
pB5 <- ggplot()+lims(x=c(-1,1))+theme_void()
pB5 <- pB5+annotate("point",x=-1,y=y,shape=16,color=c("grey90",rev(cols)),size=2)+
  annotate("text", x=-0.9,y=y,label=c("Others",rev(as.character(merged_18s.f$bacteria[match(cols,merged_18s.f$cols)]))), hjust=0)

pdf("Figure3_B.pdf",height = 3.5,width = 16)
pB1+pB2+pB3+pB4+pB5+plot_layout(nrow=1)
dev.off()
print("Completed figure3 B")

##Figure 3C
RNA_ARGs_reads <- as.character(args[14])
DNA_ARGs_reads <- as.character(args[15])

DNA_args <- read.csv(DNA_ARGs_reads,header = T,sep = '\t',row.names = 1)
RNA_args <- read.csv(RNA_ARGs_reads,header = T,sep = '\t',row.names = 1)

DNA_args1 <- DNA_args
DNA_args1 <- DNA_args1[-which(rowSums(DNA_args1)==0),]
DNA_args1[DNA_args1>1] = 1
RNA_args1 <- RNA_args
RNA_args1[RNA_args1>1] = 1


#ARM family
aro_index <- as.character(args[16])
ARO_index <- read.csv(aro_index,header = T,sep = '\t')

#ARG occurrence
overlap_args <- intersect(row.names(DNA_args1),row.names(RNA_args1))
ARGs_occurrence <- DNA_args1[overlap_args,]
for(i in overlap_args){
  for(j in 1:ncol(DNA_args1)){
    if(DNA_args1[i,j]>0 & RNA_args1[i,j]>0){ARGs_occurrence[i,j]="Both"}
    if(DNA_args1[i,j]>0 & RNA_args1[i,j]==0){ARGs_occurrence[i,j]="OnlymNGS"}
    if(DNA_args1[i,j]==0 & RNA_args1[i,j]>0){ARGs_occurrence[i,j]="OnlymtNGS"}
    if(DNA_args1[i,j]==0 & RNA_args1[i,j]==0){ARGs_occurrence[i,j]=""}
  }
}

filter_zero <- function(data){
  data <- apply(data, 1, function(x){
    length(x[x==""])/length(x)
  })
  return(data)
}
ARGs_occurrence1 <- ARGs_occurrence[filter_zero(ARGs_occurrence)<0.5,]
ARGs_occurrence1 <- ARGs_occurrence1[,c(10:62,1:9)]

#oncoPring
library(ComplexHeatmap)
mat = as.matrix(ARGs_occurrence1)
ARG_family = ARO_index$AMR.Gene.Family[match(row.names(ARGs_occurrence1),ARO_index$Model.Name)]
left_ha = rowAnnotation(GF=ARG_family,
                        col=list(GF=c("APH(3')"="#A6CEE3","APH(6)"="#A6CEE3",
                                      "rifamycin-resistant beta-subunit of RNA polymerase (rpoB)"= "#A8D1D1",
                                      "resistance-nodulation-cell division (RND) antibiotic efflux pump"="#AAD4BF",
                                      "Erm 23S ribosomal RNA methyltransferase"= "#AFDB9B",
                                      "major facilitator superfamily (MFS) antibiotic efflux pump"="#B2DF8A",
                                      "aminocoumarin resistant parY"="#C0D18D",
                                      "TEM beta-lactamase"="#CFC390",
                                      "pmr phosphoethanolamine transferase"="#DDB593",
                                      "ATP-binding cassette (ABC) antibiotic efflux pump"="#ECA796",
                                      "sulfonamide resistant sul"="#FB9A99",
                                      "tetracycline-resistant ribosomal protection protein"="#ADD8AD"))
)

bottom_ha = HeatmapAnnotation(sample_type = meta_data$sample_type[c(10:62,1:9)],
                              col = list(sample_type = c("bronchoalveolar lavage fluid sample" = "#A6CEE3", 
                                                         "cerebrospinal fluid sample" = "#B2DF8A", 
                                                         "blood sample" = "#FB9A99")))

col = c("OnlymNGS" = "#1F78B4", "OnlymtNGS" = "#E31A1C", "Both" =  "#33A02C")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  OnlymNGS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["OnlymNGS"], col = NA))
  },
  # big red
  OnlymtNGS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["OnlymtNGS"], col = NA))
  },
  # big green
  Both = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Both"], col = NA))
  }
)
column_title = "OncoPrint for Occurence of ARM genes in DNA and RNA sequencing"
heatmap_legend_param = list(title = "Prevaluence", at = c("OnlymNGS", "OnlymtNGS", "Both"), 
                            labels = c("Only in mNGS", "Only in mtNGS", "Both in mNGS and mtNGS"))

pdf("Figure3_C.pdf",height = 8,width = 16)
oncoPrint(mat,alter_fun = alter_fun, col = col,column_order = names(ARGs_occurrence1),
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,left_annotation = left_ha,
          pct_side = "right", row_names_side = "left",bottom_annotation = bottom_ha,
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)
dev.off()
print("Completed figure3 C")
print("Done")


