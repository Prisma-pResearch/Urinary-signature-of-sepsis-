# This script defines functions which prepare the data for multiple downstream analysis including predictive analytics performed in Python. 

convCEL2patID <- function(data){
  #'Changes column names from CEL file number to a unique patient ID
  pat_id <- array()
  for(i in 1:length(colnames(data))){
    pat_id[i] <- substr(colnames(data)[i], 8, 12)
  }; pat_id
  # for(i in 1:length(colnames(data))){  #only if using final prediction dataset
  #   if(grepl("_", pat_id[i]) == TRUE){pat_id[i] <- substr(pat_id[i], 1, 4)}
  # }; pat_id
  for(i in 1:length(colnames(data))){  #change the indices depending on the dataset
    if((grepl("_", pat_id[i]) == TRUE) & (grepl("N", pat_id[i]) == FALSE)){
      pat_id[i]<- paste("P0", substr(pat_id[i], 1, 3), sep = "")
    }
  }; pat_id
  for(i in 1:length(colnames(data))){  #only if you have Control patients
    if(grepl("N", pat_id[i]) == TRUE){
      pat_id[i] <- substr(pat_id[i], 1, 4)
    }
  }; pat_id
  colnames(data) <- pat_id
  return(data)
} 


getGenelist <- function(data){
  #'Get the list of unique genes present in dataset.  
  temp <- strsplit(as.character(data$Gene.Symbol), ";")
  genes_unique <- unique(unlist(temp))
  return(genes_unique)
}


limmaPrepData <- function(data){
  #'Prepares a data for LIMMA analysis 
  #'
  #'@description This function transposes the dataset, then adds an outcome column 
  limma_genes <- as.data.frame(t(data))
  limma_genes$outcome <- 1  #Patient Label is needed with gene expression value
  limma_genes$outcome[grepl("N", rownames(limma_genes))] = 0
  return(limma_genes)
}


format_forIPA <- function(data){
  #'Prepares a data for Ingenuity Pathway analysis 
  #'
  #'@description This function saves gene names as a new column. 
  ipa_data <- data
  symbol <- as.data.frame(rownames(ipa_data))  #Save the Gene Names separately in another dataframe. Dataframe is required for cbind.
  rownames(ipa_data) <- seq(1, nrow(data), 1)
  ipa_data <- cbind(symbol, ipa_data)
  return(ipa_data)
}


borutaPrepData <- function(data, ids){
  boruta_genes <- data[rownames(data) %in% rownames(ids),]
  boruta_genes <- as.data.frame(t(boruta_genes))  #Patient IDs have to be in rownames for the code in Python
  boruta_genes$outcome <- 1  #Patient Label is needed with gene expression value
  boruta_genes$outcome[grepl("N", rownames(boruta_genes))] = 0
  return(boruta_genes)
}

makeScoreSheet <- function(score, label = label, outcome = outcome, method = method, set = set){
  sheet <- cbind(score, label, outcome)
  names(sheet) <- c("Score", "Label", "Outcome")
  write.csv(sheet, file = paste ("./1 Data/Sabyasachi/ScoreSheet", method, set, ".csv", sep = "_", collapse = NULL))
}

