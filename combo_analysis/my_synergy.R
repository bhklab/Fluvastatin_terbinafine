library(abind)
library(dplyr)
library(abind)
library(robustbase)
library(synergyfinder)
# library(locfdr)
library(dplyr)
library(tidyr)



# oneill.lines <- c("KPL1","OCUBM",'T47D','EFM192B', "MDAMB436","ZR751")
source("R/reformat_sensitivity_mat.R")

getSensitivityComboMatrix <- function(viabilityFiles,x1,x2){
  
  sensitivityData <- lapply(viabilityFiles, function(x, x1,x2){
    
    print(x)
    xx <- parseComboViabilityFile(x, x1,x2)
    
    return(xx)
    
  }, x1=x1,x2=x2)
  
  sensitivityData <- abind( sensitivityData, along=1)
  
  xx <- unlist(strsplit(gsub("_[1-9]$", "", rownames(sensitivityData)), split=' (?=[^ ]+$)', perl=TRUE))
 
  
  for( i in 1:dim(sensitivityData)[1]){
    
    control <- sensitivityData[i,,][dim(sensitivityData[i,,])[1],1]
    sensitivityData[i,,] <- sensitivityData[i,,]/control
    
  }
  
  

  return(sensitivityData)
}



parseComboViabilityFile <- function(fileCombo, conc1, conc2){
  # conc1 <- x1
  #  conc2 <- x2
  
  #  fileCombo <- viabilityFiles[1]
  
  conc1 <- rev(conc1)
  viability <- read.csv(fileCombo, stringsAsFactors=FALSE, sep="\t", header=FALSE)
  
  ##this function assumes in the text file produced by the machine always there is field called "Plate:" before plate names  
  indicies <- grep("Plate", viability$V1)
  plates <- viability$V2[indicies]
  
  data.blocks <- table(is.na(viability[indicies[1] + 3, ]))["FALSE"] %/% 10#ncol(viability) %/% 10
  
  raw.sensitivity.rows <- NULL
  for(plate in plates) {
  
    xx <- unlist(strsplit(gsub("[1-9]$", "", plate), split=' (?=[^ ]+$)', perl=TRUE))
    # xx <- unlist(strsplit(gsub("[1-9]$", "", xx[1]), split=' (?=[^ ]+$)', perl=TRUE))
    cell <- paste(xx[1:(length(xx)-1)], sep=" ")
    
    drugs <- unlist(strsplit(xx[length(xx)], split='/'))
    
    if(data.blocks > 1){
      raw.sensitivity.rows <- c(raw.sensitivity.rows, do.call(c, lapply(1:data.blocks, function(x){paste(paste(cell, drugs, sep="_"), x, sep="_")})))
    }else{
      raw.sensitivity.rows <- c(raw.sensitivity.rows,cell)
    }
  }
  
  if(length(raw.sensitivity.rows)!=length(unique(raw.sensitivity.rows))){
    stop("Cell lines names are not unique")
  }
  
  
  raw.sensitivity <- array(NA, dim=c(length(raw.sensitivity.rows),length(conc1), length(conc2)), dimnames=list(raw.sensitivity.rows, conc1, conc2))
  
  for(index in indicies) {
    start <- 3
    
    plate <- viability$V2[index]
    
    xx <- unlist(strsplit(gsub("_[1-9]$", "", plate), split=' (?=[^ ]+$)', perl=TRUE))
    #    xx <- unlist(strsplit(gsub("_[1-9]$", "", xx[1]), split=' (?=[^ ]+$)', perl=TRUE))
    cell <- paste(xx[1:(length(xx)-1)], sep=" ")
    drugs <- paste(xx[(length(xx))], sep=" ")
    
    
    
    for(i in 1:length(conc1)){
      raw.sensitivity[cell,i,] <- as.numeric(do.call(c,viability[(index + i + 2), (start + 1):(start + 10)]))
      if(!is.na(viability[(index + i + 2),start])){
        raw.sensitivity[cell,i,] <- raw.sensitivity[cell,i,] - viability[(index + i + 2),start]
      }else{
        raw.sensitivity[cell,i,] <- raw.sensitivity[cell,i,] - 0.045
      }
    }
    
    
  }
  
  return(raw.sensitivity)
  
}



# Set the drug combination name
Combo <- "FLUVA_TBF"
# Define concentration levels for the two drugs:
# TBF (row drug) and Fluvastatin (col drug) as per the original experimental design
x1 <- c(0, 0.39, 1.56, 6.25, 25, 100) # concentrations for DrugR (TBF)
x2 <- c(0.0, 0.78, 1.56, 3.125, 6.25, 12.5, 25, 50, 100, 200) # concentrations for DrugC (Fluva)


files.ls <- dir("../combo_analysis/data/rawData_combo/FluvaTBF")
files.ls <- files.ls[grep("txt", files.ls)]  # select only .txt files

files.ls <- paste("../combo_analysis/data/rawData_combo/FluvaTBF/", files.ls, sep = "")

viability <- read.csv(files.ls[1], stringsAsFactors=FALSE, sep="\t", header=FALSE)
idxs <- grep("Plate", viability$V1)
plates <- viability$V2[idxs]




sensitivityData <- getSensitivityComboMatrix(viabilityFiles = files.ls, x1 = x1, x2 = x2)

sensitivityData <- sensitivityData[order(rownames(sensitivityData)), , ]

# The next line selects a subset of the matrix if needed.
# Here, we remove some concentrations or reorder the data for analysis
# In this case, we take all columns and rows 1:8 from the third dimension,
# and drop the first column if needed.
sensitivityData <- sensitivityData[ , -1, 1:8]



DrugR <- strsplit(Combo, split = "_")[[1]][2]  # TBF is DrugR (row drug)
DrugC <- strsplit(Combo, split = "_")[[1]][1]  # FLUVA is DrugC (column drug)



outfile_stats <- c()
raw_matrix <- c()
bliss_matrix <- c()


for (i in seq_len(dim(sensitivityData)[1])) {
  cellID <- rownames(sensitivityData)[i]
  mat <- sensitivityData[i,,]*100
  mat <- mat[order(as.numeric(rownames(mat))),]
  combo_row_concentration <- rownames(mat)
  combo_col_concentration <- colnames(mat)
  meta <- data.frame(drug.col = DrugC, drug.row = DrugR, concUnit = "µM", blockIDs = 1)
  data <- list(dose.response.mats = list("1" = mat), drug.pairs = meta)
  
   data_reformated = reformat_dose_response(data)
   data_reformated = data_reformated %>% select("BlockID",
   	"Response","DrugRow","DrugCol","ConcRow","ConcCol", "ConcRowUnit","ConcColUnit")

   colnames(data_reformated) <-c("BlockID","Response","DrugRow",
   	"DrugCol","ConcRow","ConcCol", "conc_r_unit","conc_c_unit")

    # data_reformated <- data_reformated %>% 
    #     filter(ConcRow<row.max) %>% 
    #     filter(ConcRow>=row.min) %>%
    #     filter(ConcCol<col.max) %>%
    #     filter(ConcCol>=col.min)

  
   
   data_reshaped = ReshapeData(data_reformated, data_type = 'viability')
   synergies <- CalculateSynergy(data = data_reshaped, correct_baseline = "non")
   zip <-  round(median(synergies$synergy_scores$ZIP_synergy, na.rm=TRUE), 4)
   
   hsa <-  round(median(synergies$synergy_scores$HSA_synergy, na.rm=TRUE), 4)
   
   bliss <-  round(median(synergies$synergy_scores$Bliss_synergy, na.rm=TRUE), 4)
   
   loewe <-  round(median(synergies$synergy_scores$ZIP_synergy, na.rm=TRUE), 4)
   
   
   result <- c(cellID, DrugR, DrugC, zip,hsa, bliss, loewe)
   
   outfile_stats <- rbind(outfile_stats, result)
 
}

print(outfile_stats)


# Format the synergy summary matrix
colnames(outfile_stats) <- c("idSample", "idDrugA", "idDrugB", "Zip","HSA","Bliss","Loewe")
print(outfile_stats)
mat_total <- as.matrix(outfile_stats[, c("Zip","HSA","Bliss","Loewe"), drop=FALSE])
class(mat_total) <- "numeric"
rownames(mat_total) <- outfile_stats[,"idSample"]
mat_total <- as.data.frame(mat_total, stringsAsFactors=FALSE)
write.csv(mat_total,"summary_synergy_trimmed.csv")