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
library(dplyr)
library(tidyr)
library(tibble)
# RNAseq <- readRDS("data/CL_data/cell_lines_expression_matrix.rds")
# print(RNAseq)

gene_mappings <- readRDS("data/CL_data/genes_ids_mappings.rds")
# print(gene_mappings)
# print(head(gene_mappings))


RNAseq <- read.csv("./data/PDO_data/gene_tpm_normalized_matrix.csv",row.names=1)

RNAseq <- RNAseq %>%t()


synergy <- read.csv("./results/organoid_synergy_scores.csv",row.names=1)
synergy <- synergy %>%t()


synergy <- synergy["Bliss",] %>% as.matrix() %>%t
rownames(synergy)<- "FLV_TBF"


gene_mappings <- readRDS("data/CL_data/genes_ids_mappings.rds")
commonSamples <- intersect(rownames(RNAseq),colnames(synergy))


gene_mappings <- gene_mappings %>% dplyr::filter(GeneBioType=='protein_coding')
# Compute correlation of each gene's expression with synergy scores
listOfAssociations <- lapply(rownames(synergy), function(y){
  geneAssociations_cor <- apply(RNAseq[commonSamples,], 2, function(x){
    
    a <- cor.test(x,synergy[y,commonSamples])
    return(c(a$estimate,a$p.value))
  })
  
  geneAssociations_cor <- t(geneAssociations_cor)
  geneAssociations_cor <- cbind(geneAssociations_cor,p.adjust(geneAssociations_cor[,2],method = "fdr"))
  colnames(geneAssociations_cor) <- c("estimate","pval","fdr")
  geneAssociations_cor <- geneAssociations_cor[order(geneAssociations_cor[,"fdr"]),]
  return(geneAssociations_cor)
})


names(listOfAssociations) <- rownames(synergy)

# Add gene symbols to the correlation results
listOfAssociations_final <- lapply(listOfAssociations, function(x){
  data.frame(x,
  "Symbol"=gene_mappings[rownames(x),"Symbol"],
  "EntrezGeneId" = gene_mappings[rownames(x),"EntrezGeneId"])
})

# Plot correlation results for one combination (FLUVA_TBF)
df <- listOfAssociations_final[[1]]
df$significant <- df$fdr < 0.05
df$significant[is.na(df$significant)] <- FALSE
write.csv(df,"./results/Organoid_Bliss_Correlation.csv")
top_genes <- df[order(df$pval), ][1:20, ]

# Scatter plot of gene correlation vs -log10 p-value, labeling top genes
pdf("./results/organoid_assoc.pdf",width = 9,height = 9)
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