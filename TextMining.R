# This function defines functions used for text mining. Text Mining in this context involves searching PubMed database for 
# papers relating our genes with "Sepsis", "Infection" or "Inflammation. 

textMining <- function(geneset){
  #keyterms<-c("sepsis AND", "infection AND", "inflammation AND")
  keyterms <- "sepsis AND"
  for (i in 1:length(keyterms))
    for(j in 1:length(geneset)){
      string<-paste("(", geneset[j], ")", sep="")#putting gene terms inside first brackets
      string<-paste(keyterms[i], string, sep=" ")#adding key terms to the front
      
      rSearch<-entrez_search(db="pubmed", term=string, retmax=90000) #entrez search
      write.csv(rSearch$ids, paste("./3 Results/Transcriptome Analysis Console/WithoutUTI_onlyRMA/Results/Text Mining_selected probes_without Lasso_only Sepsis/", geneset[j], keyterms[i], ".csv", sep="")) #write the pubmed ids in a csv file
    }
  return(rSearch)
}
