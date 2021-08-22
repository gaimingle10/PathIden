####Supplementary Figures#######
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(openxlsx)
library(ComplexHeatmap)


#FigureS1 

sample_type <- read.xlsx('supplementary/sample_type_statistic.xlsx',sheet = 1)
row.names(sample_type) <- sample_type[,1]
sample_type <- sample_type[,-1]

for(i in 1:nrow(sample_type)){
  for(j in c(1:4)){
    if(sample_type[i,j]==1){
      sample_type[i,j]=sample_type$sample_type[i]
    }else{
      sample_type[i,j]=""
    }
  }
}

hospital <- sample_type$hospital
sample_type <- sample_type[,-c(5,6)]
sample_type.t <- data.frame(t(sample_type))

col = c("cerebrospinal fluid sample" = "#377EB8",  
        "bronchoalveolar lavage fluid sample" = "#E41A1C", 
        "blood sample" =  "#4DAF4A")
top_ann = HeatmapAnnotation(Hospital = hospital,
                            col = list(Hospital = c("GuangZhou" = "#A6CEE3", 
                                                    "Peking University Third Hospital" = "#B2DF8A", 
                                                    "First Medical Center of Chinese PLA General Hospital " = "#FB9A99")))

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  `cerebrospinal fluid sample` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["cerebrospinal fluid sample"], col = NA))
  },
  # big red
  `bronchoalveolar lavage fluid sample` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["bronchoalveolar lavage fluid sample"], col = NA))
  },
  # big green
  `blood sample` = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["blood sample"], col = NA))
  }
)
column_title = "Diagram of Sequencing methods in Differential Samples"
heatmap_legend_param = list(title = "Type of Sample", 
                            at = c("cerebrospinal fluid sample", "bronchoalveolar lavage fluid sample", "blood sample"), 
                            labels = c("cerebrospinal fluid sample", "bronchoalveolar lavage fluid sample", "blood sample")
)

oncoPrint(sample_type.t,alter_fun = alter_fun, col = col,top_annotation = top_ann,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,column_order = names(sample_type.t),
          pct_side = "right", row_names_side = "left",show_column_names = TRUE,
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)

##Figure S2
#set theme for ggplot2

my_theme <- theme(
  axis.ticks.length = unit(0.4,"lines"), 
  axis.ticks = element_line(color='black'),
  axis.line = element_line(colour = "black"), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title = element_text(colour = 'black',size=8,face="bold"),
  axis.text.y=element_text(colour='black',size=8,face = "bold"),
  axis.text.x=element_text(colour='black',size=8,face = "bold"),
  plot.title=element_text(hjust = 0.5)
)

raw_data <- read.csv('supplementary/rawdata_compare.csv')
raw_data.melt <- melt(raw_data)
raw_data.melt$variable <- factor(raw_data.melt$variable,levels = c("mtNGS","mNGS"))

ggplot(raw_data.melt,aes(variable,log10(value)))+
  geom_boxplot(aes(color=variable),width=0.5,outlier.color = NA)+
  geom_jitter(aes(color=variable),position = position_jitter(0.1))+
  labs(x="",y="Log10(Raw Data Gb)",title = "Comparison of Raw Data between mtNGS and mNGS")+
  facet_wrap(.~sample_type,scales = 'free_x')+
  theme_bw()+my_theme+theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))

#FigureS3 

#FigureS4
TCID50 <- read.csv('supplementary/TCID50.csv')

ggplot(data=TCID50,aes(x=TCID50,y=Mapped_ratio))+
  labs(x="TCID50 Value",y="Mapped Ratio")+
  geom_smooth(stat="smooth",method = "lm")+
  #stat_cor(aes(label = paste(..rr.label..,..p.label..,sep = "~`,`~")),label.y = 110,color="blue")+
  stat_cor(method = "spearman",digits = 2,color="blue")+
  geom_point(stat='identity',color="#1F78B4") +
  theme_bw()+my_theme

#FigureS5
data <- read.csv('supplementary/Ucalled_readfish.csv')
data.melt <- melt(data)

ggplot(data.melt,aes(sample,value,fill=variable))+
  geom_bar(position = 'stack',stat = 'identity')+
  facet_wrap(.~Sequecing,scales = 'free_x')+
  labs(x='',y='Fraction')+scale_y_continuous(expand = c(0,0))+
  scale_fill_hue('Composition',labels=c('Host Ratio',c('Microbiota Ratio')))+
  theme_bw()+my_theme+
  theme(axis.text.x = element_text(angle = 60,hjust = 0.5,vjust = 0.5),
        legend.position = 'top')





