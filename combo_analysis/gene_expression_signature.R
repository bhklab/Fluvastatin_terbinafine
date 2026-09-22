library(survival)
library(survminer)

library(dplyr)
library(tidyr)
library(tibble)
library(survival)
data <- read.csv("./data/LegSig/organoid_cell_line_same_sign.csv")
sample_info <- read.csv("./data/TCGA_BRCA/data_clinical_sample.txt", sep="\t",skip=4)
sample_info <- sample_info %>% select("PATIENT_ID","SAMPLE_ID")

# print(sample_info[1:5,])
patient_info <- read.csv("./data/TCGA_BRCA/data_clinical_patient.txt", sep="\t",skip=4)




keep_subtypes <- c('BRCA_Her2','BRCA_LumA','BRCA_LumB','BRCA_Basal')
keep_subtypes <- c('BRCA_LumA','BRCA_LumB')
keep_subtypes <- c('BRCA_Her2') #plausible!
# keep_subtypes <- c('BRCA_Basal') 
patient_info <- patient_info %>% 
        filter(SUBTYPE %in% keep_subtypes) %>%
        # select("PATIENT_ID","OS_STATUS","OS_MONTHS","DSS_STATUS","DSS_MONTHS","DFS_STATUS","DFS_MONTHS","PFS_STATUS","PFS_MONTHS")
        select("PATIENT_ID","PFS_STATUS","PFS_MONTHS")

patient.data <- merge(patient_info,sample_info, by ="PATIENT_ID")

patient.samples <-gsub("-",".",patient.data$SAMPLE_ID)

# print(patient.samples[])



expression <- read.csv("./data/TCGA_BRCA/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",sep="\t")
expression <- expression[complete.cases(expression),] %>% filter(Hugo_Symbol!="")


# expression <- expression %>% group_by(Hugo_Symbol) %>% summarise(across(where(is.numeric), median, .names = "{col}")
# rownames(expression) <- expression$Hugo_Symbol

# stop()

cut.off <- 0.2
organoid.associations <- read.csv("./results/Organoid_Bliss_Correlation.csv")
organoid.associations <- organoid.associations %>% select(Symbol,estimate)
colnames(organoid.associations) <- c("Symbol","Organoid.Pearson")
organoid.associations <- organoid.associations[complete.cases(organoid.associations),]

cell.line.associations  <- read.csv("./results/CCL_Bliss_Correlation.csv")
cell.line.associations <- cell.line.associations %>% filter(significant==TRUE) %>% select(Symbol,estimate)
colnames(cell.line.associations) <- c("Symbol","Cell.Pearson")
cell.line.associations <- cell.line.associations[complete.cases(cell.line.associations),]
# print(organoid.associations)

expression <- expression[!duplicated(expression$Hugo_Symbol),] %>%as.data.frame()
rownames(expression)<-expression$Hugo_Symbol
expression <- expression%>%select(-Hugo_Symbol)




weight_data <- merge(cell.line.associations,organoid.associations,by='Symbol')
weight_data <- weight_data[!duplicated(weight_data$Symbol),]



weight_data$SameSign <- sign(weight_data$Cell.Pearson)==sign(weight_data$Organoid.Pearson)
weight_data <- weight_data %>% filter(Cell.Pearson>=cut.off) %>% filter(Organoid.Pearson>=cut.off)

write.csv(weight_data,"./results/merged_corr.csv")
weight_data <- weight_data %>% filter(SameSign==TRUE)


# _ <- organoid.associations %>%fil
weight_data$weights <- abs(weight_data$Cell.Pearson) / sum(abs(weight_data$Cell.Pearson))
expression <- expression %>% select(-Entrez_Gene_Id)
expression <- expression[,patient.samples]
common.genes <- intersect(rownames(expression), weight_data$Symbol)

expr_sub <- expression[common.genes, ]
common.genes <- rownames(expr_sub)
# write.csv(expr_sub[,1:5],"tumpf.csv")
weights_sub <- weight_data %>%
  filter(Symbol %in% common.genes) %>%
  arrange(match(Symbol, common.genes))
# weights_sub <- weights_sub[!duplicated(weights_sub$Symbol)]
# print(expression[1:5,1:5])
# write.csv("timpf.csv",weights_sub)

# print(weights_sub)

# print(dim(weights_sub))
# print(dim(expr_sub))

weight.vector <- as.matrix(weights_sub$weights,nrow=1,ncol=dim(weights_sub)[1])
# print(weight.vector)

# expr_sub <- t(expr_sub)


# 5. Compute scores per sample (matrix multiplication)
scores <-t(expr_sub) %*% matrix(weights_sub$weights)
scores <- as.data.frame(scores)
colnames(scores)<-"Score"
scores <- tibble::rownames_to_column(scores,"SAMPLE_ID")
# print(cl)
# col
# print(scores[1:5,])
# print(class(scores))
scores$SAMPLE_ID <-gsub("[.]","-",scores$SAMPLE_ID)

patient.data <- merge(scores, patient.data, by="SAMPLE_ID") %>%drop_na()
patient.data$PFS_STATUS <- as.numeric(sapply(strsplit(patient.data$PFS_STATUS,":"), `[`, 1))
patient.data$PFS_MONTHS <- as.numeric(patient.data$PFS_MONTHS)
# print(status)
median.score <-median(patient.data$Score)
upperQ <- quantile(patient.data$Score,2/3)
lowerQ <- quantile(patient.data$Score,1/3)

patient.data$Group<- ifelse(patient.data$Score >= median.score, "High", "Low")
patient.data$Qcut <- ifelse(patient.data$Score >=upperQ,"High",ifelse(patient.data$Score<=lowerQ,"Low","Middle"))

print(head(patient.data))
res.cox <- coxph(Surv(PFS_MONTHS, PFS_STATUS) ~ Score, data = patient.data)
print(res.cox)

km_fit <- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ Group, data = patient.data)


pdf("./results/Median_Survival_Signature.pdf")
ggsurvplot(
  km_fit,                   # The survfit object with your data
  data = patient.data,              # The dataframe used to fit the curves
  pval = TRUE,              # Displays the Log-rank test p-value
  conf.int = TRUE,          # Displays 95% confidence intervals
  risk.table = TRUE,        # Adds a "Number at Risk" table below the plot
  legend.labs = c("High", "Low"), # Customizes group names in legend
  palette = c("#E7B800", "#2E9FDF"), # Customizes curve colors
  xlab = "Time (Days)",     # Customizes X-axis label
  ggtheme = theme_minimal() # Clean background theme
)
dev.off()


patient.data <- patient.data%>% filter(Qcut %in% c("High","Low"))
km_fit <- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ Qcut, data = patient.data)


pdf("./results/Tertile_Survival_Signature.pdf")
ggsurvplot(
  km_fit,                   # The survfit object with your data
  data = patient.data,              # The dataframe used to fit the curves
  pval = TRUE,              # Displays the Log-rank test p-value
  conf.int = TRUE,          # Displays 95% confidence intervals
  risk.table = TRUE,        # Adds a "Number at Risk" table below the plot
  legend.labs = c("High", "Low"), # Customizes group names in legend
  palette = c("#E7B800", "#2E9FDF"), # Customizes curve colors
  xlab = "Time (Days)",     # Customizes X-axis label
  ggtheme = theme_minimal() # Clean background theme
)
dev.off()