#These scripts are used to set different labels for the data.

label_sepsis <- function(data){
  col <- as.data.frame(colnames(data))
  col$label <- rep(1,ncol(data))  #setting all labels to 1 initially
  rownames(col) <- col$`colnames(data)`
  col$`colnames(data)` <- NULL
  col_ind <- grep('N', colnames(data), value=TRUE)  #selecting those indices which have "N" and setting them to 0
  col[col_ind, ] <- 0
  label <- col$label; label  #Store the label separately
  return(label)
}

label_cluster <- function(data, cluster_info){
  clus_info_req <- subset(clus_info, rownames(clus_info) %in% colnames(data))  #subset out only required patients from the 'clus_info' variable
  data_total <- rbind(data, t(clus_info_req))
  label <- as.numeric(data_total[nrow(data_total), ])
  label <- replace(label, label==2, 0)  
  return(label)
}

label_CCI <- function(data, cci_outcome){
  cci_outcome_reqd <- subset(cci_outcome, rownames(cci_outcome) %in% names(data))
  data_total <- rbind(data, t(cci_outcome_reqd))
  label <- as.numeric(data_total[nrow(data_total), ])
  return(label)
}

label_AKI <- function(data, aki_outcome){
  aki_outcome_reqd <- subset(aki_outcome, rownames(aki_outcome) %in% names(data))
  names(aki_outcome_reqd) <- "outcome"
  ind <- rownames(aki_outcome_reqd)[which(aki_outcome_reqd$outcome == 5)]  #finding the name of the patient who has AKI level 5
  data[ind] <- NULL  #dropping the patient who has AKI level of 5
  aki_outcome_reqd <- subset(aki_outcome_reqd, rownames(aki_outcome_reqd) %in% names(data))  #dropping patient who has AKI
  data_total <- rbind(data, t(aki_outcome_reqd))
  label <- as.numeric(data_total[nrow(data_total), ])
  return(label)
}
