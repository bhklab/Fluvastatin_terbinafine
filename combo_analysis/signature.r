library(dplyr)
library(tibble)
library(tidyr)

data <- read.csv("./data/LegSig/organoid_cell_line_same_sign.csv")
sample_info <- read.csv("./data/TCGA_BRCA/data_clinical_sample.txt", sep="\t",skip=4)
sample_info <- sample_info %>% select("PATIENT_ID","SAMPLE_ID")

patient_info <- read.csv("./data/TCGA_BRCA/data_clinical_patient.txt", sep="\t",skip=4)
keep_subtypes <- c('BRCA_Her2','BRCA_LumA','BRCA_LumB','BRCA_Basal')
keep_subtypes <- c('BRCA_LumA','BRCA_LumB')
expression_data <- read.csv("./data/TCGA_BRCA/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",sep="\t")


patient_info <- patient_info %>% 
        filter(SUBTYPE %in% keep_subtypes) %>%
        select("PATIENT_ID","OS_STATUS","OS_MONTHS","DSS_STATUS","DSS_MONTHS","DFS_STATUS","DFS_MONTHS","PFS_STATUS","PFS_MONTHS")


print(patient_info[1:5,1:5])
print(patient_info[1])




cut_off <- 0.2

data <- data %>% 
    filter(abs(Organoid.Pearson)>=cut_off) %>%
    filter(abs(Cell.Line.Pearson)>=cut_off)

data$weights <- abs(data$Cell.Line.Pearson)/sum(abs(data$Cell.Line.Pearson))


# 1. Read your expression data (genes × samples)
expr <- read.csv("gene_tpm_organoid.csv")   # or the DF you showed
# Make Hugo_Symbol the rownames
rownames(expr) <- expr$Hugo_Symbol
expr <- expr[,-1]

# 2. Read the Bliss correlation weights
sig.weights <- read.csv("Bliss_corr.csv", row.names = 1)

# 3. Normalize weights (absolute correlation, sum = 1)
sig.weights$Weights <- abs(sig.weights$estimate) / sum(abs(sig.weights$estimate))

# 4. Align genes between expr and sig.weights
common.genes <- intersect(rownames(expr), sig.weights$Symbol)

expr_sub <- expr[common.genes, ]
weights_sub <- sig.weights %>%
  filter(Symbol %in% common.genes) %>%
  arrange(match(Symbol, common.genes))

# 5. Compute scores per sample (matrix multiplication)
scores <- as.numeric(t(expr_sub) %*% weights_sub$Weights)

# 6. Put into dataframe
score_df <- data.frame(Sample = colnames(expr_sub),
                       Score = scores)

# 7. Optional: categorize samples into High/Low
median_score <- median(score_df$Score)
score_df$Category <- ifelse(score_df$Score >= median_score, "High", "Low")

#print(score_df)

# read model list
model_list <- read.csv("model_list_all.csv")

# merge by matching Sample = modelID
merged_df <- model_list %>%
  left_join(score_df, by = c("modelID" = "Sample"))

# write out merged file
write.csv(merged_df, "model_list_test_results.csv", row.names = FALSE)

print(head(merged_df))