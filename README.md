# Urinary-signature-of-sepsis-
This repository contains full code of the paper: Discovery and validation of urinary molecular signature of early sepsis 

Paper Abstract: 
Objective: Identify alterations in gene expression unique to systemic and kidney-specific pathophysiological processes using whole-genome analyses of ribonucleic acid (RNA) isolated from the urinary cells of sepsis patients. 
Design: Prospective cohort study
Setting: Quaternary care academic hospital
Patients: 266 sepsis and 82 control patients enrolled between January 2015 and February 2018.
Interventions: Whole-genome transcriptomic analysis of mRNA isolated from the urinary cells of sepsis patients within 12 hours of sepsis onset and from control subjects. 
Measurements and Main Results: The differentially expressed probes which map to known genes were subjected to feature selection using multiple machine learning techniques to find the best subset of probes that differentiates sepsis from control subjects. Using differential expression augmented with machine learning ensembles, we identified a set of 239 genes in urine which show excellent effectiveness in classifying septic patients from those with chronic systemic disease in both internal and independent external validation cohorts. Functional analysis indexes disrupted biological pathways in early sepsis and reveals key molecular networks driving its pathogenesis.
Conclusions: We identified unique urinary gene expression profile in early sepsis. Future studies need to confirm whether this approach can complement blood transcriptomic approaches for sepsis diagnosis and prognostication.   

Brief Description of R codes. 
ModuleAnalysis.R - contains functions related to analysis of immune and kidney cell types present in urinary cells of sepsis patients. Figure 3 of paper. 
TextMining.R - contains functions for mining PubMed for certain keywords (for eg.  "sepsis", "infection" and "inflammation") and certain gene names which are found to be important in the analysis. Figure 2C of paper. 
DataPreparation.R - contains functions for manipulating dataframes in R to make them amenable for downstream machine learning modeling using python.
setLabels.R - contains functions for setting different labels to available data (for eg. sepsis/control, AKI/noAKI, CKD/noCKD, etc). 

Brief Description of Python codes. 
feature_selectores.py - contains functions related to different machine learning models (Boruta, Random Forest, SVM, Logistic regression) on the differentially expressed genes list. Each function performs Grid Search and gives the selected feature list as output. The default parameters are ones used to find the gene lists mentioned in the paper. 
preprocessing.py - contains two functions for splitting a dataframe into X and y and for splitting a gene list into upregulated and downregulated lists respectively. 
voting.py - contains one function which gives consensus voting gene lists and takes the four gene lists from the output of feature_selectors.py script. 
