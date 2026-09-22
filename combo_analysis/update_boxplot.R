require(PharmacoGx)
require(magicaxis)
library(abind)
library(robustbase)
library(Biobase)
library(synergyfinder)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(reshape)
library(snowfall)
library(GSA)
library(piano)
library(scales)
library(ggrepel)
library(dplyr)

source("./helpers/fixCellLinesNames.R")



synergies <- read.csv("./results/median_synergy_scores.csv",row.names=1)
subtypes <- read.csv("./results/ccl_subtypes.csv")
subtypes$sample <- fixCellLinesNames(subtypes$sample,"./data/CL_data/cell_annotation_all.csv")

# print(subtypes)
# print(synergies)

colnames(synergies)<- c('sample','Bliss','ZIP')
# print(synergies)

data <- merge(subtypes, synergies,by='sample') %>%  select(sample, SCMOD2, Bliss, ZIP)

scmod2_colors <- c("Basal"="#4daf4a","Her2"="#377eb8","Luminal"="#984ea3")

my_comparisons = list( c("Basal", "Luminal"), c("Luminal", "Her2"), c("Basal", "Her2") )
Subtype_col=c("Basal"="#4daf4a","Her2"="#377eb8","Luminal"="#984ea3")
#"LumA"="#e78ac3")  


##############
### Fig S6B ##
##############
# Boxplot of synergy scores grouped by SCMOD2 subtypes
pdf("./results/FigS6B_SCMOD2.pdf",width = 11,height = 6)
ggplot(data,mapping = aes(x=SCMOD2,y=Bliss,fill=SCMOD2)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_hline(yintercept = 0,lty=2,col="red") +
  geom_jitter(width = 0.2) +
  ggtitle("Boxplot of Bliss Synergy in Cell Lines") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(4.25, 5.20, 6.15))+
  stat_compare_means(label.y = 7.4) +
  scale_fill_manual(values=Subtype_col) + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x =element_text(size=14,face="bold"))
dev.off()

pdf("./results/FigS6B_SCMOD2_ZIP.pdf",width = 11,height = 6)
ggplot(data,mapping = aes(x=SCMOD2,y=ZIP,fill=SCMOD2)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_hline(yintercept = 0,lty=2,col="red") +
  geom_jitter(width = 0.2) +
  ggtitle("Boxplot of ZIP Synergy in Cell Lines") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(4.25, 5.20, 6.15))+
  stat_compare_means(label.y = 7.4) +
  scale_fill_manual(values=Subtype_col) + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x =element_text(size=14,face="bold"))
dev.off()