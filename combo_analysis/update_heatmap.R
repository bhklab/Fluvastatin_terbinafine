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
BiocManager::install("GenomeInfoDb")
load("./data/SynergyStats_Fluva_TBF.RData")
load("./data/mono_viabilities.RData")  
source("./helpers/fixCellLinesNames.R")
# source("./R/fixCellLinesNames.R")

##############################################################
### Prepare Data for Figure 4A and Figure S6B (Subtypes) #####
##############################################################

# Load breast cancer subtype annotations for the cell lines
subtypes <- read.csv("./results/ccl_subtypes.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE)


# Summarize synergy metrics (e.g., Bliss) for each combination
# Use median across replicates to get a single synergy value per cell line
listOfCombos_stat_summarized = list()

stat <- read.csv("./results/synergy_summaries.csv",row.names=1)
sampleIDs <- rownames(stat)
sampleIDs <- unique(unlist(lapply(strsplit(sampleIDs,".",fixed =TRUE),"[[",1)))


stat_summarized <- matrix(nrow = length(sampleIDs),ncol = ncol(stat),
                          dimnames = list(sampleIDs,colnames(stat)))


for (sample in sampleIDs) {
    print(sample)
    ibx <- grep(sample,rownames(stat))
    # Use median across replicates

    stat_summarized[sample,] <- colMedians(as.matrix(stat[ibx,,drop=FALSE]),na.rm = TRUE)
    if(startsWith(sample,"X")){
      print(sample)
      rownames(stat_summarized)[grep(sample,rownames(stat_summarized))] <- paste("MDAMB",gsub("X","",sample),sep = "")
    }
  }
# print(stat_summarized)
final <- fixCellLinesNames(rownames(stat_summarized),"./data/CL_data/cell_annotation_all.csv")
stn <- fixCellLinesNames(rownames(subtypes),"./data/CL_data/cell_annotation_all.csv")
# print(stn)
rownames(subtypes)<-stn
# print(final)
# stop()
rownames(stat_summarized) <- final
# print(stat_summarized)
# stop()
combo <- names(listOfCombos_stat)[1]
listOfCombos_stat_summarized[[combo]] <- stat_summarized
# print(listOfCombos_stat_summarized)

commonCells <- Reduce(intersect,lapply(listOfCombos_stat_summarized, rownames))


BlissMat_summarized <- do.call(rbind,lapply(listOfCombos_stat_summarized, function(x){
  return(x[commonCells,"Bliss"])
}))

# print(BlissMat_summarized)
# stop()

order <- names(sort(BlissMat_summarized["FLUVA_TBF",]))

# print(order)
data_heat = t(-BlissMat_summarized[,order]) 

# Extract AAC from mono treatments of Fluva and TBF
Fluva<- mono_Fluva[,"AAC"]


TBF <- mono_TBF[,"AAC"]

scmod2_colors <- c("Basal"="#4daf4a","Her2"="#377eb8","Luminal"="#e78ac3")



ref.model <- "HSA"
syn_summary <- data.frame()

ref.models <- c("ZIP","Bliss")
for(ref.model in ref.models){

  temp <- do.call(rbind,lapply(listOfCombos_stat_summarized, function(x){
  return(x[commonCells,ref.model])}))
  if(ref.model == "Bliss"){
    order <- names(sort(temp["FLUVA_TBF",]))
  }
  syn_summary <- rbind(syn_summary,temp)  
}

# print(syn_summary)


Fluva.sub <- Fluva[colnames(syn_summary)]

TBF.sub <- TBF[colnames(syn_summary)]
# order <- names(sort(syn_summary[,"Bliss"]))  # print(syn_summary)

# print(subtypes)
# print(syn_summary)
subtypes_final <- subtypes[colnames(syn_summary), ]



column_ha <- HeatmapAnnotation(
  # IntClust = subtypes_final[order, "intClust_short"],
  # PAM50 = subtypes_final[order, "PAM50"],
  SCMOD2 = subtypes_final[order,"SCMOD2"],
  FLUVA = Fluva.sub[order],
  TBF = TBF.sub[order],
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



data_heat = as.matrix(syn_summary[,order])

rownames(data_heat) <- ref.models

min_score <- floor(min(data_heat))
max_score <- ceiling(max(data_heat))


pdf("./results/Fig4_A.pdf",width = 15,height = 3)
Heatmap(data_heat,cluster_columns = FALSE,
    cluster_rows = FALSE,top_annotation = column_ha,
     heatmap_legend_param = list(
          title = "Synergy Score", at = c(min_score, max_score))
        )
dev.off()
# print("plotted")

synergy_data = as.data.frame(t(data_heat)) %>% tibble::rownames_to_column("idSample")
# print(synergy_data)
write.csv(synergy_data,"./results/median_synergy_scores.csv")




# Write out synergy scores to a CSV
data <- data.frame(t(BlissMat_summarized[,order,drop=FALSE]),check.names = FALSE)
write.table(data,"./results/synergy_scores.csv",sep = ",",row.names = TRUE,col.names = NA)

data$subtype <- subtypes[rownames(data),"SCMOD2"]
data$AAC <- Fluva[rownames(data)]

# print(data)
