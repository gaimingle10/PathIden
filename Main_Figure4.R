### Main Figure4 #####

#load Main_Figure1.R Main_Figure2.R
source("Main_Figure1.R")
source("Main_Figure2.R")

##Figrue 4A
RNA_18s <- read.csv('data/mtNGS/RNA_18s_species.txt',header = T,sep = '\t',row.names = 1)
DNA_18s <- read.csv('data/mNGS/DNA_18s_species.txt',header = T,sep = '\t',row.names = 1)

DNA_18s_data <- extract_data(DNA_18s)
RNA_18s_data <- extract_data(RNA_18s)

venn_18s <- venn_compare(RNA_18s,DNA_18s,RNA_18s_data,DNA_18s_data)

##Figure 4B
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

pB <- pB1+pB2+pB3+pB4+pB5+plot_layout(nrow=1)
pB

##Figure 4C
DNA_args <- read.csv('data/CARD/DNA_ARGs_reads.txt',header = T,sep = '\t',row.names = 1)
RNA_args <- read.csv('data/CARD/RNA_ARGs_reads.txt',header = T,sep = '\t',row.names = 1)

DNA_args1 <- DNA_args
DNA_args1 <- DNA_args1[-which(rowSums(DNA_args1)==0),]
DNA_args1[DNA_args1>1] = 1
RNA_args1 <- RNA_args
RNA_args1[RNA_args1>1] = 1


#ARM family
ARO_index <- read.csv('data/CARD/aro_index.tsv',header = T,sep = '\t')

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

oncoPrint(mat,alter_fun = alter_fun, col = col,column_order = names(ARGs_occurrence1),
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,left_annotation = left_ha,
          pct_side = "right", row_names_side = "left",bottom_annotation = bottom_ha,
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)

