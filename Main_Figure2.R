### Main Figure2 #####

#load Main_Figure1.R
source("Main_Figure1.R")

#Figure 2A
RNA_virus <- read.csv('data/mtNGS/RNA_virus_species.txt',header = T,sep = '\t',row.names = 1)
DNA_virus <- read.csv('data/mNGS/DNA_virus_species.txt',header = T,sep = '\t',row.names = 1)

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

pA <- pA1+pA2+pA3+pA4+pA5+plot_layout(nrow=1)
pA


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
#venn16s/venn23s/vennbacteria

##Figure 2B
venn_virus <- venn_compare(RNA_virus,DNA_virus,RNA_virus_data,DNA_virus_data)
venn_virus

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


ggplot(virus_data.melt,aes(bacteria,log10(value+1e-8),color=sequencing))+
  geom_boxplot(position = "dodge2")+geom_jitter(position = position_jitterdodge(0.2))+
  scale_color_manual(values = myfill[c(1,5)])+
  facet_wrap(.~sample_type,scales = 'free_x')+
  coord_flip()+labs(x='',y="Log10(Relative Abundance)")+
  theme_bw()+my_theme

