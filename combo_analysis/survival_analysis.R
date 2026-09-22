
library(survival)
library(ggplot2)
library(tibble)
library(ggsurvfit)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(survival)
library(ggplot2)
library(tibble)
patient <- read.csv("../data/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt",sep="\t",skip=4)
patient <- patient %>%
            select("PATIENT_ID","SUBTYPE","DFS_STATUS", "DFS_MONTHS","PFS_STATUS","PFS_MONTHS")


patient$PFS_STATUS <- as.numeric(substr(patient$PFS_STATUS,1,1))
patient$DFS_STATUS<- as.numeric(substr(patient$DFS_STATUS,1,1))
# patient <- patient %>%filter(SUBTYPE %in% c('BRCA_Normal'))

sample <- read.csv("../data/brca_tcga_pan_can_atlas_2018/data_clinical_sample.txt",sep="\t",skip=4)

sample <- sample %>% select("PATIENT_ID","SAMPLE_ID")

mutations <- read.csv("../data/brca_tcga_pan_can_atlas_2018/data_mutations.txt",sep="\t")
mutations$vaf <- mutations$t_alt_count/(mutations$t_ref_count + mutations$t_alt_count)
mutations <- mutations %>% filter(vaf>=0.2) %>% select("Hugo_Symbol","Tumor_Sample_Barcode")


# MUC4
# PABPC1
# MUC16
mut.vals <- mutations %>%filter(Hugo_Symbol%in%c('DMD'))

mut.samples <- unique(mut.vals$Tumor_Sample_Barcode)

print("KRING")
print(mut.samples)
# print("KRUMP")
# sample$Tumor_Sample_Barcode <- as.character(sample$Tumor_Sample_Barcode)
# f
# print(sample$Tumor_Sample_Barcode in mut.samples)
sample$mut_status <- ifelse(sample$SAMPLE_ID %in% mut.samples, 1, 0)

data <- merge(sample, patient, by='PATIENT_ID')

s1 <- survfit(Surv(PFS_MONTHS,PFS_STATUS)~mut_status,data=data)
print(summary(s1))
print(head(data))
print(s1)
print(str(s1))

png("surv.png")
ggsurvplot(
  s1,
  data = data,
  pval = TRUE,             # Adds the Log-Rank p-value to the plot
  conf.int = FALSE,      # Adds 95% confidence intervals
  risk.table = TRUE,       # Adds a risk table below the plot
  ggtheme = theme_minimal(), # Clean, modern theme
  palette = c("#E7B800", "#2E5B88") # Custom colors for the two groups
)
dev.off()
stop()


# clin.data <- read.csv("../data/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt",sep="\t",h)
weights <- read.csv("../data/organoid_cell_line_same_sign.csv")
weights <- read.csv("../data/adj_pval_ccl_organoids.csv")

weights <- weights %>%drop_na()
cutoff <- 0.2
# print(dim(weights))
# weights<- weights %>%filter(abs(Cell.Line.Pearson)>=0.2)
# print(dim(weights))
# weights <- 
expression<- read.csv(
    "../data/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt", 
    sep = "\t")

expression <- expression %>% drop_na()
print(dim(expression))

order <- sort(intersect(expression$Hugo_Symbol, weights$Symbol))
print(dim(weights))
weights <- weights %>% filter(Symbol %in% order)
print(dim(weights))
print(dim(weights))
weights <- 

