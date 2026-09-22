#####
#
# Utility functions to help read in the raw combination data
#
#
#
#  Authors: 
#     Wail Baalawi (ba3lwi) - original writing
#     James Bannon (jbannon) - renaming of variables and commenting
###


###
# Function to read in viability files from a list

getSensitivityComboMatrix <- function(viabilityFiles,x1,x2){
  
  sensitivityData <- lapply(viabilityFiles, function(x, x1,x2){
    
    
    xx <- parseComboViabilityFile(x, x1,x2)
 
    return(xx)
    
  }, x1=x1,x2=x2)

  sensitivityData <- abind( sensitivityData, along=1)
  
  xx <- unlist(strsplit(gsub("_[1-9]$", "", rownames(sensitivityData)), split=' (?=[^ ]+$)', perl=TRUE))
 
  
  for( i in 1:dim(sensitivityData)[1]){
    
    control <- sensitivityData[i,,][dim(sensitivityData[i,,])[1],1]
    sensitivityData[i,,] <- sensitivityData[i,,]/control
    
  }
  
  

  return(sensitivityData)
}






## `parse a single
parseComboViabilityFile <- function(fileCombo, conc1, conc2){
  # conc1 <- x1
  #  conc2 <- x2
  
  #  fileCombo <- viabilityFiles[1]
  # print(fileCombo)
  conc1 <- rev(conc1)
  viability <- read.csv(fileCombo, stringsAsFactors=FALSE, sep="\t", header=FALSE)
  
  ##this function assumes in the text file produced by the machine always there is field called "Plate:" before plate names  
  indicies <- grep("Plate", viability$V1)
  # print(indicies)
  plates <- viability$V2[indicies]
  # print(viability)
  data.blocks <- table(is.na(viability[indicies[1] + 3, ]))["FALSE"] %/% 10#ncol(viability) %/% 10

  raw.sensitivity.rows <- NULL
  for(plate in plates) {
  
    xx <- unlist(strsplit(gsub("[1-9]$", "", plate), split=' (?=[^ ]+$)', perl=TRUE))

    # xx <- unlist(strsplit(gsub("[1-9]$", "", x[x1]), split=' (?=[^ ]+$)', perl=TRUE))
    cell <- paste(xx[1:(length(xx)-1)], sep=" ")
    
    drugs <- unlist(strsplit(xx[length(xx)], split='/'))
    # print(data.blocks)
    # stop()
    if(data.blocks > 1){
      raw.sensitivity.rows <- c(raw.sensitivity.rows, do.call(c, lapply(1:data.blocks, function(x){paste(paste(cell, drugs, sep="_"), x, sep="_")})))
    }else{
      raw.sensitivity.rows <- c(raw.sensitivity.rows,cell)
    }
  }

  if(length(raw.sensitivity.rows)!=length(unique(raw.sensitivity.rows))){
    stop("Cell lines names are not unique")
  }
  
  
  raw.sensitivity <- array(NA, dim=c(length(raw.sensitivity.rows),length(conc1), length(conc2)), dimnames=list(raw.sensitivity.rows, conc1, conc2))
 
  # stop()
  for(index in indicies) {
    start <- 3
    
    plate <- viability$V2[index]
    
    xx <- unlist(strsplit(gsub("_[1-9]$", "", plate), split=' (?=[^ ]+$)', perl=TRUE))
    #    xx <- unlist(strsplit(gsub("_[1-9]$", "", xx[1]), split=' (?=[^ ]+$)', perl=TRUE))
    cell <- paste(xx[1:(length(xx)-1)], sep=" ")
    drugs <- paste(xx[(length(xx))], sep=" ")
    
    
    
    for(i in 1:length(conc1)){
      raw.sensitivity[cell,i,] <- as.numeric(do.call(c,viability[(index + i + 2), (start + 1):(start + 10)]))
      if(!is.na(viability[(index + i + 2),start])){
        raw.sensitivity[cell,i,] <- raw.sensitivity[cell,i,] - viability[(index + i + 2),start]
      }else{
        raw.sensitivity[cell,i,] <- raw.sensitivity[cell,i,] - 0.045 
      }
    }
    
    
  }
  
  return(raw.sensitivity)
  
}


