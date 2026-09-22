reformat_dose_response <- function(data) {
  # Extract the dose-response matrix and drug pairs
  dose_response_mat <- data$dose.response.mats$`1`
  drug_info <- data$drug.pairs
  
  # Get row and column concentrations
  conc_row <- as.numeric(rownames(dose_response_mat))
  conc_col <- as.numeric(colnames(dose_response_mat))
  
  # Initialize an empty list to store the results
  results <- list()
  
  # Iterate over the matrix and reshape into the desired format
  for (row in seq_along(conc_row)) {
    for (col in seq_along(conc_col)) {
      results <- rbind(results, data.frame(
        BlockID = 1,  # Assuming block ID is 1
        Col = col,
        Row = row,
        Response = dose_response_mat[row, col],
        Replicate = 1,  # Assuming replicate is 1
        DrugRow = drug_info$drug.row,
        DrugCol = drug_info$drug.col,
        ConcRow = conc_row[row],
        ConcCol = conc_col[col],
        ConcRowUnit = drug_info$concUnit,
        ConcColUnit = drug_info$concUnit
      ))
    }
  }
  
  # Convert results list to a data frame
  results <- as.data.frame(results)
  
  return(results)
}
