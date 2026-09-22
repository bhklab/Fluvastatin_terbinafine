# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("synergyfinder")

library(readr)
library(synergyfinder)
library(reticulate)
library(dplyr)
library(tidyr)
library(dplyr)
library(PharmacoGx)

dr1 <- 'FLV'
dr2 <- 'TBF'


# py_syn <- function(drug1, drug2 , cellline ) {
#   code <- sprintf("temp = show_report('datasets/FLV_TBF_data.csv', '%s', '%s', %d ,0, save = True, plots = True)", drug1, drug2, cellline)
#   temp <- py_run_string(code, convert = TRUE)
#   result <- py$temp[[1]]
# }

# move_res <- function(path_root,drug1, drug2){
#   today <- Sys.Date()
#   dirf_res <- paste(dirf, today, sep="/")
#   #dir.create(dirf_res, recursive = TRUE)
#   drugs <- sprintf("%s_%s", dr1, dr2)
#   dirf_res <- paste(dirf_res, drugs, sep="/")
#   dir.create(dirf_res, recursive = TRUE)
#   setwd(dirf_res)
# }

# move_root <-function(path_root){
  # setwd(dirf)
# }

# # dirf = "~/Documents/GitHub/pdo_synergy"
# move_root(dirf)

HTdataset <- read_csv("./data/PDO_data/pdo_synergy/datasets/FLV_TBF_data.csv")
colnames(HTdataset) <- c("date", "conc1", "conc2", "response","repeats", "drug1", "drug2", "conc_unit1", "conc_unit2", "block_id")

print(HTdataset)

aacs <-c()
blck.ids <- c()
reps <- c()
drugs <- c()
for(bid in unique(HTdataset$block_id)){
  
  for(drug in c('FLV','TBF')){
  
  
  if(drug=='FLV'){
    # looking at conc2 (TBF) to 0
    temp <- HTdataset %>% filter(block_id==bid,conc2 ==0)
  } else {
    #looking at FLV
    temp <- HTdataset %>% filter(block_id==bid,conc1 ==0)
  }
   
  
  
  

  
  for(rr in unique(temp$repeats)){
    if(drug=='FLV'){
      dr <- temp %>% filter(repeats==rr) %>% arrange(conc1)
       if (dim(dr)[1]>=3){
      
      aac <- 1-computeAUC(dr$conc1,dr$response)/100
      aacs <- c(aacs,aac)
      blck.ids <- c(blck.ids,bid)
      reps <- c(reps, rr)
      drugs <- c(drugs, drug)
    }
    }
    else{
      dr <- temp %>% filter(repeats==rr) %>% arrange(conc2)
       if (dim(dr)[1]>=3){
      
      aac <- 1-computeAUC(dr$conc2,dr$response)/100
      aacs <- c(aacs,aac)
      blck.ids <- c(blck.ids,bid)
      reps <- c(reps, rr)
      drugs <- c(drugs, drug)
    }
    }
  }
}
}

aac.results <- data.frame(block=blck.ids,replicate = reps, aac = aacs, drug = drugs)
write.csv(aac.results,file="./results/pdo_sensitivity.csv")
# move_res(dirf, dr1, dr2)

pairs <- subset(HTdataset,  drug1==dr1 & drug2==dr2 & block_id != 64)
single <- subset(HTdataset,  drug1==dr1 & conc2==0 & block_id != 64)
single$drug2 = dr2
SynData <- rbind(pairs, single)

res <- ReshapeData(
  data = SynData,
  data_type = "viability",
  impute = TRUE,
  impute_method = NULL,
  iteration = 10,
  seed = 1)

cellLines <-res$drug_pairs$block_id


res <- CalculateSynergy(
  data = res,
  method = c("ZIP", "Bliss"),
  Emin = 0,
  Emax = 1,
  iteration =10,
  correct_baseline = "non")



res <- CalculateSensitivity(
  data = res,
  correct_baseline = "non", 
  iteration = 10
)

# PlotDoseResponse(
#   data = res,
#   summary_statistic = "mean",
#   block_ids = cellLines,
#   drugs = c(1,2),
#   statistic  = "ci",
#   save_file = TRUE,
#   file_type = "png"
# )




# PlotSynergy(
#   data = res,
#   type = "heatmap",
#   method = "ZIP",
#   block_ids = cellLines,
#   drugs = c(1,2),
#   save_file = TRUE,
#   file_type = "png", 
#   grid = FALSE,
#   heatmap_text_label_size_scale = 0, 
#   high_value_color = "#D55E00",
#   low_value_color = "#0072B2",
# )







report <- res$drug_pairs
print("Starting PDO writes")
write.csv(res$drug_pairs, file="./results/pdo_pairs.csv" ,row.names=FALSE, quote=FALSE) 
# write.csv(res$synergy_scores, file="./results/pdo_synergy.csv", dr1, dr2, row.names=FALSE, quote=FALSE) 
print("PDO Synergy Done")