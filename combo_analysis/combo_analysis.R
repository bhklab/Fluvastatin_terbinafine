#############################
### Required R Packages #####
#############################
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

#############################
### Important: Directories ###
#############################

# Set the working directory to the "FluvaTBF_repo-main" folder.
# Adjust the path below as needed for your local environment.
setwd("Fluvastatin_terbinafine/combo_analysis/")

# Set whether you want to recalculate the synergy values from raw data.
# If TRUE, synergy values are computed from scratch using "calculateSynergy_combo.R".
# If FALSE, precomputed data is loaded from "data/SynergyStats_Fluva_TBF.RData".
PerformSynergyCalcFromRaw = TRUE

########################################
### Calculating/Loading Synergy Data ###
########################################

if(PerformSynergyCalcFromRaw){
  # Source script to calculate synergy from raw dose-response data
  source("R/calculateSynergy_combo.R")
}else{
  # Load precomputed synergy statistics
  load("data/SynergyStats_Fluva_TBF.RData",verbose = TRUE)
}

# Source script to standardize cell line names across experiments
source("R/fixCellLinesNames.R")

############################################################
### Compute Fluvastatin Monotherapy Effects from Combo Data #
############################################################

# For each combination experiment, extract the Fluvastatin-only response values
# and calculate IC50 and AAC (Area Above Curve).

listOfCombos_MonoFluva <- lapply(names(listOfCombos), function(combo){
  
  # Extract mono treatment results by taking the first column of the dose-response matrix (Fluva axis)
  monoResults <- do.call(rbind,lapply(names(listOfCombos[[combo]]$Bliss),function(sample){
    unlist(listOfCombos[[combo]]$Bliss[[sample]]$dose.response.mats[[1]][,1])
  }))
  
  rownames(monoResults) <- names(listOfCombos[[combo]]$Bliss)
  colnames(monoResults) <- rownames(listOfCombos[[combo]]$Bliss[[1]]$dose.response.mats[[1]])
  
  # Average across replicates using median
  uniq_samples <- unique(unlist(lapply(strsplit(rownames(monoResults)," "),"[[",1)))
  monoResults_final <- do.call(rbind,lapply(uniq_samples, function(x){
    ibx <- grep(x,x = rownames(monoResults))
    return(colMedians(monoResults[ibx,,drop=FALSE]))
  }))
  
  rownames(monoResults_final) <- uniq_samples
  # Fix cell line names
  rownames(monoResults_final) <- fixCellLinesNames(rownames(monoResults_final),"data/CL_data/cell_annotation_all.csv")
  
  # Compute IC50 and AAC from the concentration-response data
  fluvaMono <- lapply(1:dim(monoResults_final)[1],function(x){
    ic50 <- computeIC50(concentration = as.numeric(colnames(monoResults_final)[-1]),
                        viability = monoResults_final[x,-1],
                        viability_as_pct = TRUE)
    aac <- computeAUC(concentration = as.numeric(colnames(monoResults_final)[-1]),
                      viability = monoResults_final[x,-1],
                      viability_as_pct = TRUE)/100
    return(c(ic50,aac))
  })
  
  fluvaMono <- do.call(rbind,fluvaMono)
  rownames(fluvaMono) <- rownames(monoResults_final)
  colnames(fluvaMono) <- c("IC50","AAC")
  
  return(fluvaMono)
})

names(listOfCombos_MonoFluva) <- names(listOfCombos)

#############################################################
### Compute TBF (Tobrafenib) Monotherapy Effects from Combo #
#############################################################

listOfCombos_MonoTBF <- lapply(names(listOfCombos), function(combo){
  
  # Extract TBF-only response values by taking the first row of the dose-response matrix
  monoResults <- do.call(rbind,lapply(names(listOfCombos[[combo]]$Bliss),function(sample){
    unlist(listOfCombos[[combo]]$Bliss[[sample]]$dose.response.mats[[1]][1,])
  }))
  
  rownames(monoResults) <- names(listOfCombos[[combo]]$Bliss)
  colnames(monoResults) <- colnames(listOfCombos[[combo]]$Bliss[[1]]$dose.response.mats[[1]])
  
  # Average across replicates using median
  uniq_samples <- unique(unlist(lapply(strsplit(rownames(monoResults)," "),"[[",1)))
  monoResults_final <- do.call(rbind,lapply(uniq_samples, function(x){
    ibx <- grep(x,x = rownames(monoResults))
    return(colMedians(monoResults[ibx,,drop=FALSE]))
  }))
  
  rownames(monoResults_final) <- uniq_samples
  # Fix cell line names
  rownames(monoResults_final) <- fixCellLinesNames(rownames(monoResults_final),"data/CL_data/cell_annotation_all.csv")
  
  # Compute IC50 and AAC for TBF mono data
  TBFMono <- lapply(1:dim(monoResults_final)[1],function(x){
    ic50 <- computeIC50(concentration = as.numeric(colnames(monoResults_final)[-1]),
                        viability = monoResults_final[x,-1],
                        viability_as_pct = TRUE)
    aac <- computeAUC(concentration = as.numeric(colnames(monoResults_final)[-1]),
                      viability = monoResults_final[x,-1],
                      viability_as_pct = TRUE)/100
    return(c(ic50,aac))
  })
  
  TBFMono <- do.call(rbind,TBFMono)
  rownames(TBFMono) <- rownames(monoResults_final)
  colnames(TBFMono) <- c("IC50","AAC")
  
  return(TBFMono)
})

names(listOfCombos_MonoTBF) <- names(listOfCombos)

##############################################################
### Prepare Data for Figure 4A and Figure S6B (Subtypes) #####
##############################################################

# Load breast cancer subtype annotations for the cell lines
subtypes <- read.csv("data/CL_data/ cell_subtypes.tsv",header = TRUE,row.names = 1,stringsAsFactors = FALSE,sep = "\t")
subtypes <- subtypes[,c("PAM50","SCMOD2","intClust")]
subtypes$intClust_short = sapply(strsplit(subtypes$intClust, " "), `[`, 1)

# Summarize synergy metrics (e.g., Bliss) for each combination
# Use median across replicates to get a single synergy value per cell line
listOfCombos_stat_summarized = list()

for (i in 1:length(listOfCombos_stat)) {
  combo <- names(listOfCombos_stat)[i]
  stat <- listOfCombos_stat[[combo]]
  
  sampleIDs <- rownames(stat)
  sampleIDs <- unique(unlist(lapply(strsplit(sampleIDs,".",fixed =TRUE),"[[",1)))
  
  stat_summarized <- matrix(nrow = length(sampleIDs),ncol = ncol(stat),
                            dimnames = list(sampleIDs,colnames(stat)))
  
  for (sample in sampleIDs) {
    ibx <- grep(sample,rownames(stat))
    # Use median across replicates
    stat_summarized[sample,] <- colMedians(as.matrix(stat[ibx,,drop=FALSE]),na.rm = TRUE)
    
    # Standardize cell line names
    if(startsWith(sample,"X")){
      rownames(stat_summarized)[grep(sample,rownames(stat_summarized))] <- paste("MDAMB",gsub("X","",sample),sep = "")
    }
  }
  
  final <- fixCellLinesNames(rownames(stat_summarized),"data/CL_data/cell_annotation_all.csv")
  rownames(stat_summarized) <- final
  listOfCombos_stat_summarized[[combo]] <- stat_summarized
}

# Identify common cell lines and combine synergy data across combinations
commonCells <- Reduce(intersect,lapply(listOfCombos_stat_summarized, rownames))
BlissMat_summarized <- do.call(rbind,lapply(listOfCombos_stat_summarized, function(x){
  return(x[commonCells,"Bliss"])
}))

# Extract AAC from mono treatments of Fluva and TBF
Fluva <- listOfCombos_MonoFluva[[1]][,"AAC"]
TBF <- listOfCombos_MonoTBF[[1]][,"AAC"]

# Subset data and order columns by FLUVA_TBF synergy
Fluva.sub <- Fluva[colnames(BlissMat_summarized)]
TBF.sub <- TBF[colnames(BlissMat_summarized)]
order <- names(sort(BlissMat_summarized["FLUVA_TBF",]))

subtypes_final <- subtypes[colnames(BlissMat_summarized), ]

# Define subtype color schemes
pam50_colors <- c("Basal"="#4daf4a","Her2"="#377eb8","LumB"="#984ea3","LumA"="#e78ac3","Normal"="yellow2")  
scmod2_colors <- c("Basal"="#4daf4a","HER2"="#377eb8","LumB"="#984ea3","LumA"="#e78ac3")
intclust_types <- unique(na.omit(subtypes_final$intClust_short))
intclust_colors = setNames(topo.colors(length(intclust_types)), intclust_types)

# Create column annotations for the heatmap
# This includes subtype annotations and mono-treatment sensitivities
column_ha <- HeatmapAnnotation(
  IntClust = subtypes_final[order, "intClust_short"],
  PAM50 = subtypes_final[order, "PAM50"],
  SCMOD2 = subtypes_final[order, "SCMOD2"],
  FLUVA = Fluva.sub[order],
  TBF = TBF.sub[order],
  annotation_legend_param = list(
    FLUVA = list(
      title = "FLUVA",
      at = c(0, 0.5, 1),
      labels = c("Insensitive", "", "Sensitive")
    ),
    TBF = list(
      title = "TBF",
      at = c(0, 0.5, 1),
      labels = c("Insensitive", "", "Sensitive")
    )
  ),
  col = list(
    IntClust = intclust_colors,
    PAM50 = pam50_colors,
    SCMOD2 = scmod2_colors,
    FLUVA = colorRamp2(c(0,0.5,1), (rev(c("red","white","blue")))),
    TBF = colorRamp2(c(0,0.5,1), (rev(c("red","white","blue"))))
  )
)

##############
### Fig 4A ###
##############
# Plot a heatmap of Bliss synergy values across cell lines for FLUVA_TBF

data_heat = t(-BlissMat_summarized[,order])  # Negative sign to invert synergy sign if needed
rownames(data_heat) <- "Bliss"

pdf("./results/Fig4_A.pdf",width = 15,height = 3)
Heatmap(data_heat,cluster_columns = FALSE,cluster_rows = FALSE,top_annotation = column_ha,
        col = colorRamp2(c(-3,0,3),(c("aquamarine3","white","wheat3"))),
        heatmap_legend_param = list(
          title = "Bliss", at = c(-3, 3), 
          labels = c("Synergy", "Antagonism")
        ))
dev.off()

# Write out synergy scores to a CSV
data <- data.frame(t(BlissMat_summarized[,order,drop=FALSE]),check.names = FALSE)
write.table(data,"./results/synergy_scores.csv",sep = ",",row.names = TRUE,col.names = NA)

data$subtype <- subtypes[rownames(data),"SCMOD2"]
data$AAC <- Fluva[rownames(data)]

# Preparing data for boxplots (Fig S6B)
BlissMat_summarized_all_final <- data.frame(t(BlissMat_summarized),
                                            "intClust"=subtypes[colnames(BlissMat_summarized),"intClust_short"],
                                            stringsAsFactors = FALSE)
BlissMat_summarized_all_final$cellID <- rownames(BlissMat_summarized_all_final)
BlissMat_summarized_all_final <- melt(BlissMat_summarized_all_final, id=c("cellID","intClust")) 
colnames(BlissMat_summarized_all_final) <- c("cellID", "intClust", "DrugCombo", "Synergy")
BlissMat_summarized_all_final$intClust <- as.factor(BlissMat_summarized_all_final$intClust)

Subtype_col=intclust_colors

##############
### Fig S6B ##
##############
# Boxplot of synergy scores grouped by IntClust subtype
my_comparisons = list( c("iC1", "iC10"), c("iC1", "iC4"))

pdf("./results/FigS6B_intClust.pdf",width = 11,height = 6)
ggplot(BlissMat_summarized_all_final,mapping = aes(x=intClust,y=Synergy,fill=intClust)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_hline(yintercept = 0,lty=2,col="red") +
  geom_jitter(width = 0.2) +
  ggtitle("FLUVA_TBF") +
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 7.4) +
  scale_fill_manual(values=Subtype_col) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x =element_text(size=14,face="bold"))
dev.off()

# Repeat similar analysis for SCMOD2 subtypes
BlissMat_summarized_all_final <- data.frame(t(BlissMat_summarized),
                                            "SCMOD2"=subtypes[colnames(BlissMat_summarized),"SCMOD2"],
                                            stringsAsFactors = FALSE)
BlissMat_summarized_all_final$cellID <- rownames(BlissMat_summarized_all_final)
BlissMat_summarized_all_final <- melt(BlissMat_summarized_all_final, id=c("cellID","SCMOD2")) 
colnames(BlissMat_summarized_all_final) <- c("cellID", "SCMOD2", "DrugCombo", "Synergy")
BlissMat_summarized_all_final$SCMOD2 <- as.factor(BlissMat_summarized_all_final$SCMOD2)

my_comparisons = list( c("Basal", "LumB"), c("LumB", "HER2"), c("Basal", "HER2") )
Subtype_col=c("Basal"="#4daf4a","HER2"="#377eb8","LumB"="#984ea3","LumA"="#e78ac3")  

##############
### Fig S6B ##
##############
# Boxplot of synergy scores grouped by SCMOD2 subtypes
pdf("./results/FigS6B_SCMOD2.pdf",width = 11,height = 6)
ggplot(BlissMat_summarized_all_final,mapping = aes(x=SCMOD2,y=Synergy,fill=SCMOD2)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_hline(yintercept = 0,lty=2,col="red") +
  geom_jitter(width = 0.2) +
  ggtitle("FLUVA_TBF") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(4.25, 5.20, 6.15))+
  stat_compare_means(label.y = 7.4) +
  scale_fill_manual(values=Subtype_col) + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x =element_text(size=14,face="bold"))
dev.off()

#######################################
# Gene Expression Correlations with Synergy
#######################################

RNAseq <- readRDS("data/CL_data/cell_lines_expression_matrix.rds")
commonSamples <- intersect(rownames(RNAseq),colnames(BlissMat_summarized))
gene_mappings <- readRDS("data/CL_data/genes_ids_mappings.rds")

# Compute correlation of each gene's expression with synergy scores
listOfAssociations <- lapply(rownames(BlissMat_summarized), function(y){
  geneAssociations_cor <- apply(RNAseq[commonSamples,], 2, function(x){
    a <- cor.test(x,BlissMat_summarized[y,commonSamples])
    return(c(a$estimate,a$p.value))
  })
  
  geneAssociations_cor <- t(geneAssociations_cor)
  geneAssociations_cor <- cbind(geneAssociations_cor,p.adjust(geneAssociations_cor[,2],method = "fdr"))
  colnames(geneAssociations_cor) <- c("estimate","pval","fdr")
  geneAssociations_cor <- geneAssociations_cor[order(geneAssociations_cor[,"fdr"]),]
  return(geneAssociations_cor)
})

names(listOfAssociations) <- rownames(BlissMat_summarized)

# Add gene symbols to the correlation results
listOfAssociations_final <- lapply(listOfAssociations, function(x){
  data.frame(x,"Symbol"=gene_mappings[rownames(x),"Symbol"])
})

# Plot correlation results for one combination (FLUVA_TBF)
df <- listOfAssociations_final[[1]]
df$significant <- df$fdr < 0.05
df$significant[is.na(df$significant)] <- FALSE

top_genes <- df[order(df$pval), ][1:20, ]

# Scatter plot of gene correlation vs -log10 p-value, labeling top genes
pdf("./results/Fig4_B.pdf",width = 9,height = 9)
ggplot(df, aes(x = estimate, y = -log10(pval), color = significant)) +
  geom_point() +
  geom_text_repel(data = top_genes, 
                  aes(label = Symbol),
                  size = 4,
                  color = "black",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps=20) +
  theme_bw() +
  labs(x = "Estimate", 
       y = "-log10(p-value)",
       title = "Gene Expression Association with Synergy",
       color = "FDR < 0.05")
dev.off()

#######################################
# Pathway Enrichment (GSEA) on Gene Correlations
#######################################

nPerm <- 10000
gsc1 <- loadGSC("data/pathwaysInfo/h.all.v6.2.symbols.gmt")

listOfPathways_enriched <- lapply(names(listOfAssociations), function(combo){
  
  sigAssoc_final_all <- listOfAssociations[[combo]]
  genesSymbols <- gene_mappings[rownames(sigAssoc_final_all),"Symbol"]
  ibx <- which(duplicated(genesSymbols))
  
  gene_stat_levels <- sigAssoc_final_all[-ibx,"estimate"]
  names(gene_stat_levels) <- genesSymbols[-ibx]
  
  gene_stat_levels <- gene_stat_levels[!is.na(gene_stat_levels)]
  gene_stat_levels <- sort(gene_stat_levels,decreasing = TRUE)
  
  set.seed(3425)
  # Run GSEA using piano
  gseaRes <- piano::runGSA(geneLevelStats = gene_stat_levels,gsc = gsc1,adjMethod = "none",
                           nPerm = nPerm,geneSetStat = "fgsea",ncpus = 4)
  
  gseaResSummary <- piano::GSAsummaryTable(gseaRes)
  gseares1 <- gseaResSummary[,c("Name","Genes (tot)","Stat (dist.dir)",
                                "p (dist.dir.up)","p (dist.dir.dn)",
                                "Genes (up)","Genes (down)")]
  gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE)
  # Replace p=0 with p=1/(nPerm+1)
  gseares1[gseares1[,"pval"] == 0,"pval"] <- 1/(nPerm+1)
  
  gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr")
  gseares1 <- gseares1[,-(4:5)]
  gseares1 <- gseares1[order(gseares1$FDR,gseares1$pval,-abs(gseares1$`Stat (dist.dir)`),na.last = TRUE),]
  
  return(gseares1)
})

names(listOfPathways_enriched) <- names(listOfAssociations)

listOfPathways_enriched_sig <- lapply(listOfPathways_enriched, function(x){
  gseares1Sig <- x[!is.na(x$FDR) & x$FDR<=0.05,]
  gseares1Sig <- gseares1Sig[order((gseares1Sig$`Stat (dist.dir)`)),]
  return(gseares1Sig)
})

# For FLUVA_TBF combination
Pathways_enriched_sig = listOfPathways_enriched_sig[[1]]

# Plot top pathways
df = Pathways_enriched_sig
top_n <- 33
df_top <- df[order(df$pval), ][1:top_n, ]
df_top$Name <- factor(df_top$Name, levels = rev(df_top$Name))
colnames(df_top) <- c("Name","Count", "Stat", "Genes (up)", "Genes (down)","pval", "FDR")
df_top = df_top[order(df_top$Stat),]

pdf("./results/Fig4_C.pdf",height = 15,width = 8)
ggplot(df_top, aes(x = Stat, y = reorder(Name, Stat))) +
  geom_point(aes(size = Count, color = FDR)) +
  scale_color_gradient(high = "blue", low = "red") +
  theme_bw() +
  xlab("Stat") + 
  ylab("Pathway") +
  ggtitle("") +
  theme(axis.text.y = element_text(size = 10))
dev.off()


#######################################
# SQLE Expression vs Synergy/AAC Correlation
#######################################

# SQLE gene ID: ENSG00000104549.11

# SYNERGY vs SQLE
common <- intersect(names(RNAseq[,'ENSG00000104549.11']),names(BlissMat_summarized[1,]))
x <- RNAseq[common,'ENSG00000104549.11']
y <- BlissMat_summarized["FLUVA_TBF",common]

# Pearson correlation
cor_test_pearson <- cor.test(x, y, method = "pearson")
corr_value <- round(cor_test_pearson$estimate, 2)
p_value <- signif(cor_test_pearson$p.value, 3)
annot_text <- paste0("r = ", corr_value, "\np = ", p_value)

df <- data.frame(x = x, y = y)

pdf("./results/SQLE_Synergy.pdf")
ggplot(df, aes(x = x, y = y)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = Inf, y = Inf, label = annot_text, hjust = 1.1, vjust = 2, size = 5) +
  theme_bw(base_size = 14) +
  labs(title = "SQLE gene expression vs FLUVA_TBF Synergy",
       x = "log2(TPM+1)",
       y = "Bliss") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
dev.off()

# SQLE vs FLUVA AAC
common <- intersect(names(RNAseq[,'ENSG00000104549.11']),names(Fluva))
x <- RNAseq[common,'ENSG00000104549.11']
y <- Fluva[common]

cor_test_pearson <- cor.test(x, y, method = "pearson")
corr_value <- round(cor_test_pearson$estimate, 2)
p_value <- signif(cor_test_pearson$p.value, 3)
annot_text <- paste0("r = ", corr_value, "\np = ", p_value)

df <- data.frame(x = x, y = y)

pdf("./results/SQLE_FLUVA.pdf")
ggplot(df, aes(x = x, y = y)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = Inf, y = Inf, label = annot_text, hjust = 1.1, vjust = 2, size = 5) +
  theme_bw(base_size = 14) +
  labs(title = "SQLE gene expression vs FLUVA",
       x = "log2(TPM+1)",
       y = "AAC") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
dev.off()

# SQLE vs TBF AAC
common <- intersect(names(RNAseq[,'ENSG00000104549.11']),names(TBF))
x <- RNAseq[common,'ENSG00000104549.11']
y <- TBF[common]

cor_test_pearson <- cor.test(x, y, method = "pearson")
corr_value <- round(cor_test_pearson$estimate, 2)
p_value <- signif(cor_test_pearson$p.value, 3)
annot_text <- paste0("r = ", corr_value, "\np = ", p_value)

df <- data.frame(x = x, y = y)

pdf("./results/SQLE_TBF.pdf")
ggplot(df, aes(x = x, y = y)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate("text", x = Inf, y = Inf, label = annot_text, hjust = 1.1, vjust = 2, size = 5) +
  theme_bw(base_size = 14) +
  labs(title = "SQLE gene expression vs TBF",
       x = "log2(TPM+1)",
       y = "AAC") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
dev.off()

######################################################
# Proteomic Data (RPPA) Analysis for Associations #####
######################################################

RPPA_UHN <- read.csv("data/CL_data/breast_rppa_UHN_NormLog2.tsv",sep = "\t",header = FALSE,stringsAsFactors = FALSE)

# Extract normalized RPPA data matrix
RPPA_UHN_final <- (RPPA_UHN[9:dim(RPPA_UHN)[1],11:dim(RPPA_UHN)[2]])
rownames(RPPA_UHN_final) <- RPPA_UHN[9:dim(RPPA_UHN)[1],7]
colnames(RPPA_UHN_final) <- RPPA_UHN[5,11:dim(RPPA_UHN)[2]]

RPPA_UHN_final <- as.matrix(RPPA_UHN_final)
class(RPPA_UHN_final) <- "numeric"

gene_names_RPPA <- as.character(RPPA_UHN[3,11:dim(RPPA_UHN)[2]])
names(gene_names_RPPA) <- colnames(RPPA_UHN_final)

cells <- fixCellLinesNames(rownames(RPPA_UHN_final),"data/CL_data/cell_annotation_all.csv")
rownames(RPPA_UHN_final) <- cells

RPPA_UHN_mapping <- as.character(RPPA_UHN[3,11:dim(RPPA_UHN)[2]])
names(RPPA_UHN_mapping) <- RPPA_UHN[5,11:dim(RPPA_UHN)[2]]

commonSamples <- intersect(rownames(RPPA_UHN_final),colnames(BlissMat_summarized))

# Correlate RPPA protein levels with synergy
listOfAssociations_RPPA <- lapply(rownames(BlissMat_summarized), function(y){
  
  RPPA_Associations_cor <- apply(RPPA_UHN_final[commonSamples,], 2, function(x){
    a <- cor.test(x,BlissMat_summarized[y,commonSamples])
    return(c(a$estimate,a$p.value))
  })
  
  RPPA_Associations_cor <- t(RPPA_Associations_cor)
  RPPA_Associations_cor <- cbind(RPPA_Associations_cor,p.adjust(RPPA_Associations_cor[,2],method = "fdr"))
  colnames(RPPA_Associations_cor) <- c("estimate","pval","fdr")
  RPPA_Associations_cor <- RPPA_Associations_cor[order(RPPA_Associations_cor[,"fdr"]),]
})

names(listOfAssociations_RPPA) <- rownames(BlissMat_summarized)

df <- as.data.frame(listOfAssociations_RPPA[[1]])
df$significant <- df$fdr < 0.05
df$significant[is.na(df$significant)] <- FALSE
df$Symbol <- gene_names_RPPA[rownames(df)]
top_genes <- df[order(df$pval), ][1:6, ]

# Scatter plot of protein correlation vs -log10 p-value
pdf("./results/FigSupp_Proteomics.pdf",width = 9,height = 9)
ggplot(df, aes(x = estimate, y = -log10(pval), color = significant)) +
  geom_point() +
  geom_text_repel(data = top_genes, 
                  aes(label = Symbol),
                  size = 4,
                  color = "black",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps=20) +
  theme_bw() +
  labs(x = "Estimate", 
       y = "-log10(p-value)",
       title = "Protein Expression Association with Synergy",
       color = "FDR < 0.05")
dev.off()

##############################################
# Generate heatmaps of median synergy scores #
# for selected cell lines (HCC1143,HCC1937,  #
# MDAMB231,SUM149)                           #
##############################################

samples_test = c("HCC1143","HCC1937","MDAMB231","SUM149")

pdf("./results/All_Samples_Heatmaps.pdf", width = 7, height = 6)

for (sample in samples_test) {
  # Find replicates for the given sample
  ibx <- grep(sample, names(listOfCombos$FLUVA_TBF$Bliss))
  
  # Extract score matrices for all replicates
  score_list <- lapply(listOfCombos$FLUVA_TBF$Bliss[ibx], function(x) x$scores[[1]])
  
  # Combine into a 3D array: rows x columns x replicates
  arr <- array(unlist(score_list),
               dim = c(nrow(score_list[[1]]), ncol(score_list[[1]]), length(score_list)))
  
  # Median across replicates for each concentration pair
  median_mat <- apply(arr, c(1, 2), median)
  dimnames(median_mat) <- dimnames(score_list[[1]])
  
  # Modify matrix by removing first row/column and reversing column order as needed
  median_mat_tmp <- median_mat
  median_mat_tmp[1, ] <- NA
  median_mat_tmp[, 1] <- NA
  median_mat_tmp <- median_mat_tmp[, rev(colnames(median_mat_tmp))]
  
  # Compute overall median synergy value
  overall_median <- median(median_mat_tmp, na.rm = TRUE)
  
  # Color scale from blue (-30) to white (0) to red (30)
  col_fun <- colorRamp2(c(-30, 0, 30), c("blue", "white", "red"))
  
  # Heatmap for this sample
  ht <- Heatmap(t(median_mat_tmp),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                col = col_fun,
                name = "Bliss",
                column_title = "FLUVA (uM)",
                row_title = "TBF (uM)",
                column_title_side = "bottom",
                row_title_side = "left",
                row_names_side = "left",
                heatmap_legend_param = list(title = "Bliss")
  )
  
  # Add a title with sample name and overall median value
  draw(ht, column_title = paste0(sample, "\nMedian: ", round(overall_median, 2)))
}

dev.off()
