# This script defines functions which prepare the data for multiple downstream analysis including predictive analytics performed in Python. 

convCEL2patID <- function(data){
  patID <-array()
  for(i in 1:length(colnames(data))){
    patID[i] <- substr(colnames(data)[i], 8, 12)
  }; patID
  # for(i in 1:length(colnames(data))){                                             #only if using final prediction dataset
  #   if(grepl("_", patID[i]) == TRUE){patID[i] <- substr(patID[i], 1, 4)}
  # }; patID
  for(i in 1:length(colnames(data))){                                               #change the indices depending on the dataset
    if((grepl("_", patID[i]) == TRUE) & (grepl("N", patID[i]) == FALSE)){patID[i]<- paste("P0", substr(patID[i], 1, 3), sep = "")}
  }; patID
  for(i in 1:length(colnames(data))){                                             #only if you have Control patients
    if(grepl("N", patID[i]) == TRUE){patID[i] <- substr(patID[i], 1, 4)}
  }; patID
  colnames(data) <- patID
  return(data)
} 


getGenelist <- function(data){
  temp <- strsplit(as.character(data$Gene.Symbol), ";")
  genes_unique <- unique(unlist(temp))
  return(genes_unique)
}


limmaPrepData <- function(data){
  limma_genes<-as.data.frame(t(data))
  limma_genes$Outcome <- 1 #Patient Label is needed with gene expression value
  limma_genes$Outcome[grepl("N", rownames(limma_genes))] = 0
  
  return(limma_genes)
}


format_forIPA <- function(data){
  ipa_data <- data
  symbol<-as.data.frame(rownames(ipa_data)) #Save the Gene Names separately in another dataframe. Dataframe is required for cbind.
  rownames(ipa_data)<-seq(1, nrow(data), 1)
  ipa_data <- cbind(symbol, ipa_data)
  return(ipa_data)
}


borutaPrepData <- function(data, ids){
  boruta_genes <- data[rownames(data) %in% rownames(ids),]
  boruta_genes<-as.data.frame(t(boruta_genes)) #Patient IDs have to be in rownames for the code in Python
  boruta_genes$Outcome <- 1 #Patient Label is needed with gene expression value
  boruta_genes$Outcome[grepl("N", rownames(boruta_genes))] = 0
  
  return(boruta_genes)
}

makeScoreSheet <- function(score, label = label, Outcome = Outcome, method = method, set = set){
  sheet <- cbind(score, label, Outcome)
  names(sheet) <- c("Score", "Label", "Outcome")
  write.csv(sheet, file = paste ("./1 Data/Sabyasachi/ScoreSheet", method, set, ".csv", sep = "_", collapse = NULL))
}

