#####
#
#
#	Code to compute the synergy values 
#	for Fluvastatin and TBF across
#	multiple breast cancer cell lines. 
#
#
#
#	Author: James Bannon
#	Email: bhklab.jamesbannon@gmail.com
#	Github: jbannon
#
######

library(abind)
library(dplyr)
library(robustbase)
library(synergyfinder)
library(tidyr)

source("./helpers/file_utils.R")
source("./helpers/reformat_sensitivity_mat.R")
####
#
#	Raw Data
#
###

# TBF := Tobrafenib 
# FLUVA: Fluvastatin

# stop()
raw.data.path <- "./data/rawData_combo/FluvaTBF/"

# Names and screening concentrations of the drugs used in the combo
row.Drug <- "TBF"
row.Conc <- c(0, 0.39, 1.56, 6.25, 25, 100) 

col.Drug <- "FLUVA"
col.Conc <- c(0.0, 0.78, 1.56, 3.125, 6.25, 12.5, 25, 50, 100, 200) # concentrations for DrugC (Fluva)


# Collect the files with the raw 
# combination screening files
files.ls <- dir(raw.data.path)
files.ls <- files.ls[grep("txt", dir(raw.data.path))]  # select only .txt files
files.ls <- paste(raw.data.path, files.ls, sep = "")

# viability <- read.csv(files.ls[1], stringsAsFactors=FALSE, sep="\t", header=FALSE)
# idxs <- grep("Plate", viability$V1)
# plates <- viability$V2[idxs]


sensitivityData <- getSensitivityComboMatrix(
    viabilityFiles = files.ls, 
    x1 = row.Conc, x2 = col.Conc)

sensitivityData <- sensitivityData[order(rownames(sensitivityData)), , ]


sensitivityData <- sensitivityData[ , -1, 1:8]



outfile_stats <- c()
raw_matrix <- c()
zip.regions <- tibble()

for (i in seq_len(dim(sensitivityData)[1])) {
  cellID <- rownames(sensitivityData)[i]
  
  
  mat <- sensitivityData[i,,]*100

  mat <- mat[order(as.numeric(rownames(mat))),]
  combo_row_concentration <- rownames(mat)
  combo_col_concentration <- colnames(mat)
  meta <- data.frame(drug.col = col.Drug, drug.row = row.Drug, concUnit = "µM", blockIDs = 1)
  data <- list(dose.response.mats = list("1" = mat), drug.pairs = meta)
  
  data_reformated = reformat_dose_response(data)
 
  data_reformated = data_reformated %>% select("BlockID",
   	"Response","DrugRow","DrugCol","ConcRow","ConcCol", "ConcRowUnit","ConcColUnit")

  colnames(data_reformated) <-c("BlockID","Response","DrugRow",
   	"DrugCol","ConcRow","ConcCol", "conc_r_unit","conc_c_unit")

  # data_reformated <- data_reformated %>% 
  #     filter(ConcRow<25) %>% 
  #     filter(ConcCol<25)
  data_reshaped = ReshapeData(data_reformated, data_type = 'viability')
 
  
  synergies <- CalculateSynergy(data = data_reshaped,
   method = c("Bliss", "HSA", "ZIP", "Loewe"),
   correct_baseline = "non")

  
  synergies$synergy_scores <- synergies$synergy_scores %>% filter(conc1>0, conc2>0)

  zip <-  round(median(synergies$synergy_scores$ZIP_synergy, na.rm=TRUE), 4) 
  hsa <-  round(median(synergies$synergy_scores$HSA_synergy, na.rm=TRUE), 4)
  bliss <-  round(median(synergies$synergy_scores$Bliss_synergy, na.rm=TRUE), 4)
  loewe <-  round(median(synergies$synergy_scores$Loewe_synergy, na.rm=TRUE), 4)
  

  synergies$response <- synergies$response %>% filter(conc1>0, conc2>0)



  long.mat <- synergies$synergy_scores %>% 
    select(block_id, conc1, conc2, ZIP_ref, Bliss_ref, Loewe_ref, HSA_ref)%>% 
    merge(synergies$response) %>%    
    select(-response_origin)


  long.mat$idSample <- cellID

  raw_matrix <- rbind(raw_matrix,long.mat)
  positive.zip.regions <- synergies$synergy_scores %>% 
    filter(ZIP_synergy>0) %>% 
    select(block_id,conc1, conc2, ZIP_synergy, Bliss_synergy, Loewe_synergy, HSA_synergy)
  
  positive.zip.regions$idSample <- cellID
  positive.zip.regions$drug1 <- synergies$drug_pairs$drug1[1]
  positive.zip.regions$drug2 <- synergies$drug_pairs$drug2[1]
  zip.regions <- rbind(zip.regions,positive.zip.regions)

  result <- c(cellID, row.Drug, col.Drug, bliss,hsa, zip, loewe)
  
  outfile_stats <- rbind(outfile_stats, result)
 
}



# Format the synergy summary matrix
colnames(outfile_stats) <- c("idSample", "idDrugA", "idDrugB", "Bliss", "HSA", "ZIP", "Loewe")

mat_total <- as.matrix(outfile_stats[, c("Bliss", "HSA", "ZIP", "Loewe"), drop=FALSE])
class(mat_total) <- "numeric"
rownames(mat_total) <- outfile_stats[,"idSample"]
mat_total <- as.data.frame(mat_total, stringsAsFactors=FALSE)

zip.regions <- as.data.frame(zip.regions, stringsAsFactors = FALSE)

raw_matrix <- as.data.frame(raw_matrix, stringsAsFactors = FALSE)
write.csv(raw_matrix,"./results/long_synergy_scores.csv")
write.csv(mat_total, "./results/synergy_summaries.csv")
write.csv(zip.regions, "./results/zip_regions.csv")