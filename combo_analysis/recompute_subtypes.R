 # nolint # nolint: indentation_linter.
 #!/usr/bin/env Rscript

# --- install libraries ---
BiocManager::install("genefu")
BiocManager::install("org.Hs.eg.db")

# --- Load required libraries ---
suppressPackageStartupMessages({
  library(genefu)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})


write.dir <- "./results/"
organoid.meta <- c("organoid", "gene_tpm_normalized_matrix.csv", "./data/PDO_data/")

ccl.meta <- c("ccl", "cell_lines_expression_matrix.rds","./data/CL_data/")

meta <- list(organoid = organoid.meta, ccl = ccl.meta)

for(source in names(meta)){
    print(source)
    info <- meta[[source]]
    
    data.source <- info[1]
    
    base.dir <- info[3]

    input.file <- paste0(base.dir,info[2])
    print(input.file)
    
    output.file <- paste0(write.dir,data.source,"_subtypes.csv")

    if(data.source == 'organoid'){
    

        expr <- read.csv(input.file, row.names = 1, check.names = FALSE)
        expr_log2 <- log2(expr + 1)
    } else{
        expr <- readRDS(input.file)
        expr <- 2**expr - 0.01
        expr.log2.tpm.p1 <- log2(expr+1)
        
        expr_log2 <- t(expr.log2.tpm.p1)
    }



    # --- Step 2: Map Ensembl → Entrez IDs ---
    rownames(expr_log2) <- sub("\\..*", "", rownames(expr_log2))
    entrez_map <- mapIds(
    org.Hs.eg.db,
    keys       = rownames(expr_log2),
    keytype    = "ENSEMBL",
    column     = "ENTREZID",
    multiVals  = "first")
    
    # --- Step 3: Build annotation data frame ---
    annot <- data.frame(
    probe         = rownames(expr_log2),
    EntrezGene.ID = entrez_map,
    row.names     = rownames(expr_log2))

    expr_t <- t(expr_log2)


    # --- Step 5: Load SCMOD2 model ---
    if ("scmod2.robust" %in% data(package = "genefu")$results[, "Item"]) {
    message("✅ scmod2.robust object found → loading robust centroids")
    data(scmod2.robust)
    } else {
    message("⚠️ scmod2.robust not found → using classic pam50")
    data(scmod2)
    }


    sbt.model <- "scmod2"

    # --- Step 6: Run PAM50 subtyping ---
    subtype_res <- molecular.subtyping(
        sbt.model  = sbt.model,   # always "pam50"
        data       = expr_t,
        annot      = annot,
        do.mapping = TRUE
)


# print(subtype_res)
    
  #Select samples pertaining to Basal Subtype
    Basals<-names(which(subtype_res$subtype == "ER-/HER2-"))
    #Select samples pertaining to HER2 Subtype
    HER2s<-names(which(subtype_res$subtype == "HER2+"))
    #Select samples pertaining to Luminal Subtypes
    LuminalB<-names(which(subtype_res$subtype == "ER+/HER2- High Prolif"))
    LuminalA<-names(which(subtype_res$subtype == "ER+/HER2- Low Prolif"))
    
    # print(source)
    # print(LuminalB)
    # print(LuminalA)
    # print(Basals)
    
    
    
    
    # --- Step 7: Save results ---
    subtype_calls <- data.frame(
    sample  = rownames(expr_t),
    subtype = subtype_res$subtype,
    SCMOD2 = NA
    )

    
    subtype_calls[which(subtype_calls$sample %in% c(LuminalB,LuminalA)),]$SCMOD2 <- "Luminal"
    
    # subtype_calls[which(subtype_calls$sample %in% LuminalA),]$SCMOD2 <- "Luminal"
    
    subtype_calls[which(subtype_calls$sample %in% Basals),]$SCMOD2 <- "Basal"
    subtype_calls[which(subtype_calls$sample %in% HER2s),]$SCMOD2 <- "Her2"
    # # print(subtype_res)
    # stop()
    print(subtype_calls)

    write.csv(subtype_calls, output.file, row.names = FALSE)

cat("✅ Subbtyping complete. Results saved to:", output.file, "\n")

}
