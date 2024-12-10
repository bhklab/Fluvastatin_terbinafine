#############################
### Required R Packages #####
#############################
require(PharmacoGx)
require(magicaxis)
library(abind)
library(robustbase)
library(Biobase)
library(synergyfinder)

# Source custom functions needed for this analysis
source("R/synergySurf.R")
source("R/fixCellLinesNames.R")
source("R/SynVolume.r")
source("R/AnalyzeCombo_as_function.R")
source("R/reformat_sensitivity_mat.R")

#############################
### Initialize Data Lists ####
#############################

# listOfCombos will store the synergy results for each combination
# listOfCombos_stat will store summarized synergy stats (e.g., median Bliss scores per sample)
listOfCombos <- list()
listOfCombos_stat <- list()

############################################################
### Parsing Sensitivity Raw Files for FLUVA_TBF Combo ######
############################################################

# Set the drug combination name
Combo <- "FLUVA_TBF"

# Define concentration levels for the two drugs:
# TBF (row drug) and Fluvastatin (col drug) as per the original experimental design
x1 <- c(0, 0.39, 1.56, 6.25, 25, 100) # concentrations for DrugR (TBF)
x2 <- c(0.0, 0.78, 1.56, 3.125, 6.25, 12.5, 25, 50, 100, 200) # concentrations for DrugC (Fluva)

# List files containing raw viability data for Fluva+TBF combinations
files.ls <- dir("../combo_analysis/data/rawData_combo/FluvaTBF")
files.ls <- files.ls[grep("txt", files.ls)]  # select only .txt files
files.ls <- paste("../combo_analysis/data/rawData_combo/FluvaTBF/", files.ls, sep = "")

# Parse the raw files and assemble them into a 3D array:
# Dimensions: samples x DrugR_concentrations x DrugC_concentrations
sensitivityData <- getSensitivityComboMatrix(viabilityFiles = files.ls, x1 = x1, x2 = x2)

# Sort the rows by sample name
sensitivityData <- sensitivityData[order(rownames(sensitivityData)), , ]

# The next line selects a subset of the matrix if needed.
# Here, we remove some concentrations or reorder the data for analysis
# In this case, we take all columns and rows 1:8 from the third dimension,
# and drop the first column if needed.
sensitivityData <- sensitivityData[ , -1, 1:8]

############################################################
### Calculate Bliss Synergy Scores for Each Sample #########
############################################################

DrugR <- strsplit(Combo, split = "_")[[1]][2]  # TBF is DrugR (row drug)
DrugC <- strsplit(Combo, split = "_")[[1]][1]  # FLUVA is DrugC (column drug)

outfile_stats <- c()
raw_matrix <- c()
bliss_matrix <- c()
BlissList <- list()

# Loop through each sample (row) in sensitivityData
for (i in seq_len(dim(sensitivityData)[1])) {
  cellID <- rownames(sensitivityData)[i]
  
  # Extract viability matrix for the current sample and convert to percentage
  mat <- sensitivityData[i,,]*100
  # Ensure row ordering by concentration
  mat <- mat[order(as.numeric(rownames(mat))),]
  
  combo_row_concentration <- rownames(mat)
  combo_col_concentration <- colnames(mat)
  
  # Construct meta info for synergy calculation
  meta <- data.frame(drug.col = DrugC, drug.row = DrugR, concUnit = "µM", blockIDs = 1)
  data <- list(dose.response.mats = list("1" = mat), drug.pairs = meta)
  
  # Try to calculate Bliss synergy scores using synergyfinder functions
  bliss_score <- NA
  bliss_mat <- NA
  
  bliss_exe <- tryCatch({ 
    data_reformated = reformat_dose_response(data)
    data_reshaped = ReshapeData(data_reformated, data.type = 'viability')
    bliss <- CalculateSynergy(data = data_reshaped, method = "Bliss", correction = FALSE)
    bliss$dose.response.mats[[1]] <- data$dose.response.mats$`1`
    bliss
  }, error = function(err) { 
    print("Bliss error"); 
    bliss <- NA
    return(bliss)
  })
  
  # If Bliss calculation worked, summarize synergy score as median of inner matrix
  if(all(!is.na(bliss_exe))) {
    bliss_score <- round(median(bliss_exe$scores[[1]][-1,-1], na.rm=TRUE), 4)
    bliss_mat <- bliss_exe$scores[[1]]
    BlissList[[i]] <- bliss_exe
  } else {
    # If no synergy data, fill with NA
    bliss_mat <- matrix(rep(NA, length(c(mat))), ncol = ncol(mat), nrow = nrow(mat))
    rownames(bliss_mat) <- combo_row_concentration
    colnames(bliss_mat) <- combo_col_concentration
    BlissList[[i]] <- NULL
  }
  
  # Store summary stats
  result <- c(cellID, DrugR, DrugC, bliss_score)
  outfile_stats <- rbind(outfile_stats, result)
  
  # Store raw viability data in a long format
  mat_out <- data.frame(
    name = rep(cellID, length(c(mat))),
    drugA = rep(DrugR, length(c(mat))),
    drugB = rep(DrugC, length(c(mat))),
    concA = as.numeric(rep(colnames(mat), each=nrow(mat))),
    concB = as.numeric(rep(rownames(mat), ncol(mat))),
    experiments = as.numeric(c(mat))
  )
  raw_matrix <- rbind(raw_matrix, mat_out)
  
  # Store Bliss synergy data in a long format
  bliss_out <- data.frame(
    name = rep(cellID, length(c(bliss_mat))),
    drugA = rep(DrugR, length(c(bliss_mat))),
    drugB = rep(DrugC, length(c(bliss_mat))),
    concA = as.numeric(rep(c(combo_col_concentration), each=nrow(bliss_mat))),
    concB = as.numeric(rep(c(combo_row_concentration), ncol(bliss_mat))),
    bliss = as.numeric(c(bliss_mat))
  )
  bliss_matrix <- rbind(bliss_matrix, bliss_out)
}

# Assign sample names to BlissList
names(BlissList) <- rownames(sensitivityData)

# Format the synergy summary matrix
colnames(outfile_stats) <- c("idSample", "idDrugA", "idDrugB", "Bliss")

mat_total <- as.matrix(outfile_stats[, "Bliss", drop=FALSE])
class(mat_total) <- "numeric"
rownames(mat_total) <- outfile_stats[,"idSample"]
mat_total <- as.data.frame(mat_total, stringsAsFactors=FALSE)

# Store results in global lists
listOfCombos[[Combo]] <- list("Bliss"=BlissList)
listOfCombos_stat[[Combo]] <- mat_total

# Save the synergy results for the FLUVA_TBF combination
save(listOfCombos, listOfCombos_stat, file = "data/SynergyStats_Fluva_TBF.RData")
