#this is the main script. It will call functions from auxiliary scripts and do some basic operations. 

setwd(the directory where your codes are saved) #please change path to suit your code. 
source("setLabels.R")
source("ModuleAnalysis.R")
source("DataPreparation.R")
source("TextMining.R")

####### Bioconductor Libraries required 
library(Biobase)
library(limma)
library(mygene)
library(ComplexHeatmap)

####### CRAN Libraries required 
library(ggplot2)  
library(circlize)  
library(reshape2)  
library(EnvStats)  
library(devtools)  
library(rentrez)  
library(plyr)  
library(e1071)
library(tidyr)
library(dplyr)
library(plyr) 

setwd(the directory where your Data and Results are saved)

exp_matrix_all <- read.csv("./1 Data/TAC software/Expression matrices/Expression matrix_all_noUTI.csv", 
                           header = TRUE, 
                           row.names = 1)
fc_pval_all <- read.csv("./1 Data/TAC software/Fold Change_P Val/Limma_all.csv", 
                        header = TRUE, 
                        row.names = 1)
fc_pval_known <- fc_pval_all[!(as.numeric(fc_pval_all$Gene.Symbol) == 1),]
exp_matrix_all <- convCEL2patID(exp_matrix_all)  #DataPreparation.R script. Converts column names from CEL file name to patient ID. 

label = label_sepsis(exp_matrix_all)  #setLabels.R script
outcome<-as.factor(label)
outcome<-as.character(outcome) 
for (i in 1:length(outcome)){  #Converting 0/1 to Control/Case. This is required by many functions
  if (outcome[i] == 1) {
    outcome[i]<- 'Sepsis'
  } else {
    outcome[i] <- 'Vascular'
  }
}
exp_matrix_known <- exp_matrix_all[rownames(fc_pval_known),]

#Cellular-Deconvolution. All functions here can be found in ModuleAnalysis.R script.
iris<-read.csv("./1 Data/IRIS_data.csv", header=TRUE)
iris_reduced <- as.data.frame(iris[iris$Cell.Specificity %in% 
                                   c("Neutrophil", "Dendritic Cell", "B Cell", 
                                     "T Cell", "NK Cell", "Monocyte"), 
                                   c("Name", "Cell.Specificity")])
write.csv(iris_reduced, "./1 Data/iris_genesymbol.csv")
a2s.obj_clean <- alias2symbol_transform(iris_reduced)
write.csv(a2s.obj_clean, "./1 Data/iris_genesymbol_alias2symbol.csv")
fc_pval_known_unnest <- data_unnest(fc_pval_known)
probe_cell_uniquePID <- probe2cell(fc_pval_known_unnest, 
                                   a2s.obj_clean)
write.csv(probe_cell_uniquePID, "./1 Data/Probe2Cell_avgexpression.csv")
total_iris <- analysis_ready_data(exp_matrix_known, probe_cell_uniquePID)
write.csv(total_iris, "./1 Data/IRIS_probes_allinfo.csv")
heatmap_sorted(exp_matrix_iris, total_iris)
aggdata_iris_trunc <- aggregate_data(total_iris, exp_matrix_iris)
write.csv(aggdata_iris_trunc, "./1 Data/IRIS_probe_TAC_forstatisticaltest.csv")
avg_exp_plt(aggdata_iris)
dir_reg_plt(total_iris)

#Post-feature selection tasks
exp_matrix_diff <- read.csv("./1 Data/TAC software/Expression matrices/Expression matrix_diffexp_noUTI.csv", 
                            header = TRUE, 
                            row.names = 1)
fc_pval_diff <- read.csv("./1 Data/TAC software/Fold Change_P Val/Limma_FDR0.01_logFC1.csv", 
                         header = TRUE, 
                         row.names = 1)
fc_pval_diff_known <- fc_pval_diff[!(as.numeric(fc_pval_diff$Gene.Symbol) == 1),]
voted_probes <- read.csv("./1 Data/Voted_probes.csv", header = TRUE, row.names = 1)
fc_pval_vote <- fc_pval_diff_known[voted_probes$Features, ]
write.csv(fc_pval_vote, "./1 Data/FC_pVal_votedprobes.csv")
genes_unique <- getGenelist(fc_pval_vote) #DataPreparation.R script. 

#Text Mining of these selected genes. Find function in TextMining.R 
textMining(genes_unique)

#Preparing data for Ingenuity Pathway Analysis (IPA)
newdata_ipa <- format_forIPA(fc_pval_vote)
write.table(newdata_ipa, "./1 Data/IPA/IPA_TAC_RMA_selectedprobes.txt", 
            quote = FALSE, 
            sep="\t")











