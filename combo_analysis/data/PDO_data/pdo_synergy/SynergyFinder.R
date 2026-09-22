# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("synergyfinder")

library(readr)
library(synergyfinder)
library(reticulate)

dr1 <- 'FLV'
dr2 <- 'TBF'


py_syn <- function(drug1, drug2 , cellline ) {
  code <- sprintf("temp = show_report('datasets/FLV_TBF_data.csv', '%s', '%s', %d ,0, save = True, plots = True)", drug1, drug2, cellline)
  temp <- py_run_string(code, convert = TRUE)
  result <- py$temp[[1]]
}

move_res <- function(path_root,drug1, drug2){
  today <- Sys.Date()
  dirf_res <- paste(dirf, today, sep="/")
  #dir.create(dirf_res, recursive = TRUE)
  drugs <- sprintf("%s_%s", dr1, dr2)
  dirf_res <- paste(dirf_res, drugs, sep="/")
  dir.create(dirf_res, recursive = TRUE)
  setwd(dirf_res)
}

move_root <-function(path_root){
  setwd(dirf)
}

dirf = "~/Documents/GitHub/pdo_synergy"
move_root(dirf)

HTdataset <- read_csv("./datasets/FLV_TBF_data.csv")
colnames(HTdataset) <- c("date", "conc1", "conc2", "response","repeats", "drug1", "drug2", "conc_unit1", "conc_unit2", "block_id")

move_res(dirf, dr1, dr2)

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

PlotDoseResponse(
  data = res,
  summary_statistic = "mean",
  block_ids = cellLines,
  drugs = c(1,2),
  statistic  = "ci",
  save_file = TRUE,
  file_type = "png"
)


#png(file=sprintf("%s_%s_sens-syn.png", dr1, dr2))
#PlotSensitivitySynergy(
#  data = res,
#  plot_synergy = "ZIP",
#  show_labels = TRUE,
#  dynamic = FALSE
#)
#dev.off()

PlotSynergy(
  data = res,
  type = "heatmap",
  method = "ZIP",
  block_ids = cellLines,
  drugs = c(1,2),
  save_file = TRUE,
  file_type = "png", 
  grid = FALSE,
  heatmap_text_label_size_scale = 0, 
  high_value_color = "#D55E00",
  low_value_color = "#0072B2",
)

move_root(dirf)


move_res(dirf, dr1, dr2)
report <- res$drug_pairs
write.csv(report, file=sprintf("%s_%s_report.csv", dr1, dr2), row.names=FALSE, quote=FALSE) 
write.csv(res$synergy_scores, file=sprintf("%s_%s_data.csv", dr1, dr2), row.names=FALSE, quote=FALSE) 