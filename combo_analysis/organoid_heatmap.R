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
library(tibble)
library(dplyr)
library(tidyr)

subtypes <- read.csv("./results/organoid_subtypes.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE)
synergies <- read.csv("./results/pdo_pairs.csv")
map_ids <- read.csv("./data/PDO_data/map_block_organoidID.csv")
synergies <- synergies %>% select('block_id','ZIP_synergy', 'Bliss_synergy')
synergies <- merge(synergies, map_ids,by='block_id') %>% select(-block_id) %>% filter(PDO_name!='BPTO95')
colnames(synergies)<- c("ZIP","Bliss","PDO_name")
rownames(synergies) <- synergies$PDO_name
synergies <- synergies %>% select(-PDO_name)
synergies <- t(synergies)
# synergies <- synergies[,-c('BPTO95')]
sensitivities <- read.csv("./results/pdo_sensitivity.csv",row.names=1)
scmod2_colors <- c("Basal"="#4daf4a","Her2"="#377eb8","Luminal"="#e78ac3")
print(map_ids)
print(sensitivities)
colnames(sensitivities)<-c('block_id','replicate','aac','drug')
sensitivities <- merge(sensitivities,map_ids,by='block_id') %>% select(-block_id,-replicate) %>% filter(PDO_name!='BPTO95')

FLV <- sensitivities %>% filter(drug=='FLV')
TBF <- sensitivities %>% filter(drug=='TBF')

FLV <- FLV %>% group_by(PDO_name)%>%summarise(AAC = median(aac,na.rm=TRUE)) %>% as.data.frame() 
TBF <- TBF %>% group_by(PDO_name)%>%summarise(AAC = median(aac,na.rm=TRUE))%>% as.data.frame()

rownames(FLV)<-FLV$PDO_name
# FLV <- FLV %>% select(-PDO_name) %>%t()
# rownames(FLV)<-NULL
FLV.final <- as.numeric(FLV$AAC)
TBF.final <- as.numeric(TBF$AAC)
names(TBF.final) <- TBF$PDO_name
names(FLV.final)<-FLV$PDO_name



order <- names(sort(synergies["Bliss",]))


# print(subtypes[order,"SCMOD2"])


# print(synergies)
# stop()
column_ha <- HeatmapAnnotation(
  # IntClust = subtypes_final[order, "intClust_short"],
  # PAM50 = subtypes_final[order, "PAM50"],
  SCMOD2 = subtypes[order,"SCMOD2"],
  FLUVA = FLV.final[order],
  TBF = TBF.final[order],
  annotation_legend_param = list(
    FLUVA = list(
      title = "FLUVA",
      at = c(0, 0.5, 1)
      # labels = c("Insensitive", "", "Sensitive")
    ),
    TBF = list(
      title = "TBF",
      at = c(0, 0.5, 1)
      # labels = c("Insensitive", "", "Sensitive")
    )
  ),
  col = list(
    # IntClust = intclust_colors,
    # PAM50 = pam50_colors,
    SCMOD2 = scmod2_colors,
    FLUVA = colorRamp2(c(0,0.5,1), (rev(c("red","white","blue")))),
    TBF = colorRamp2(c(0,0.5,1), (rev(c("red","white","blue"))))
  )
)





data_heat = as.matrix(synergies[,order])
# 
# rownames(data_heat) <- ref.models


min_score <- floor(min(data_heat))
max_score <- ceiling(max(data_heat))


pdf("./results/organoid_heatmap.pdf",width = 15,height = 3)
Heatmap(data_heat,cluster_columns = FALSE,
    cluster_rows = FALSE,top_annotation = column_ha,
     heatmap_legend_param = list(
          title = "Synergy Score", at = c(min_score, max_score))
        )
dev.off()


synergy_data = as.data.frame(t(data_heat)) %>% tibble::rownames_to_column("idSample")
# print(synergy_data)
write.csv(synergy_data,"./results/organoid_synergy_scores.csv",row.names = FALSE)

