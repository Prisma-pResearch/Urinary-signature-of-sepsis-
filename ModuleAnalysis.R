#this script runs functions for immune/kidney analysis

alias2symbol_transform <- functions(data){
  a2s.obj <- alias2SymbolUsingNCBI(data$Name, "./4 Resources/Homo_sapiens.gene_info", required.columns = c("GeneID","Symbol"))
  a2s.obj <- cbind(a2s.obj, data)
  a2s.obj_clean <- a2s.obj[!is.na(a2s.obj$GeneID),]
  colnames(a2s.obj_clean)[2] <- "Gene.Symbol"
  return(a2s.obj_clean)
}


data_unnest <- function(data){
  data$ProbeIds <- rownames(data)
  data_unnest <- data %>% 
    mutate(Gene.Symbol = strsplit(as.character(Gene.Symbol), ";")) %>% 
    unnest(Gene.Symbol)
  return(data_unnest)
}


probe2cell <- function(fc_pval_unnest, alia2sym){
  probe_cell <- merge(fc_pval_unnest,alia2sym,by="Gene.Symbol")
  probe_cell <- probe_cell[order(probe_cell$Cell.Specificity),]
  probe_cell_uniquePID <- probe_cell[!duplicated(probe_cell$ProbeIds),] #removing duplicated ProbeIds
  rownames(probe_cell_uniquePID) <- probe_cell_uniquePID$ProbeIds
  probe_cell_uniquePID$ProbeIds <- NULL
  return(probe_cell_uniquePID)
}


analysis_ready_data <- function(exp_mat, probe2cell){
  exp_matrix <- exp_mat[rownames(probe2cell),] #expression values of only IRIS genes
  total <- cbind(probe2cell, exp_mat)
  return(total)
}


heatmap_sorted <- function(exp_mat, total){
  data_sepsis <- exp_mat[,grepl("P", colnames(exp_mat))]; data_control <- exp_mat[,grepl("N", colnames(exp_mat))]
  data_sepsis_n<-data_sepsis[, order(colSums(-data_sepsis))]; data_control_n<-data_control[, order(colSums(-data_control))]
  index_order = c(names(data_control_n), names(data_sepsis_n))
  data_hm = exp_mat[,index_order]
  y <- c(replicate(ncol(data_control), "Vascular"), replicate(ncol(data_sepsis), "Sepsis"))
  
  mat = as.matrix(data_hm) 
  base_mean = rowMeans(mat)
  mat_scaled = t(apply(mat, 1, scale))
  type = colnames(data_hm) 
  ha = HeatmapAnnotation(data.frame(total$Cell.Specificity), 
                         col = list(total.Cell.Specificity = c("B Cell" = "mediumblue", "Dendritic Cell" = "chartreuse3", "Monocyte" = "slategray3", "Neutrophil" =  "pink2", "NK Cell" = "orange", "T Cell" = "red")),
                         annotation_legend_param = list(total.Cell.Specificity = list(title = "Immune cell type", title_gp = gpar(fontsize = 12), 
                                                                                           labels_gp = gpar(fontsize = 10))))
  setEPS() 
  postscript("Heatmap_sortedboth_immune.eps")
  Heatmap(t(mat_scaled), name = "Expression", col = colorRamp2(c(-3, 0, 3), c("green", "black", "red")),
          top_annotation = ha, top_annotation_height = unit(4, "mm"), column_names_side = "top",
          show_row_names = TRUE, show_column_names = FALSE, show_row_dend = TRUE, show_column_dend = FALSE,
          cluster_rows = FALSE, cluster_columns = FALSE)+
    Heatmap(y, name = "", width = unit(2.5, "mm"), col = c("Sepsis" = "indianred", "Vascular" = "cornflowerblue"))
  dev.off()
}


aggregate_data <- function(total, exp_mat){
  data <- cbind(total$Cell.Specificity, exp_mat); colnames(data)[1] <- "Cell.Specificity" #Add Cell.specificity to the front of expression values
  aggdata <-aggregate(data, list(total$Cell.Specificity), FUN=mean) 
  aggdata$Cell.Specificity<-aggdata$Group.1; aggdata$Group.1<-NULL
  rownames(aggdata)<-aggdata$Cell.Specificity; aggdata$Cell.Specificity<-NULL
  aggdata <- aggdata[,c(colnames(aggdata)[grepl("P", colnames(aggdata))], colnames(aggdata)[grepl("N", colnames(aggdata))])] #reordering colnames as "Cell.Specificity" -> "Sepsis" -> "Vascular"
  return(aggdata)
}


avg_exp_plt <- function(aggdata){
  bplot.obj<-as.data.frame(t(aggdata)); 
  colNew <- c(replicate(table(grepl("P", colnames(aggdata)))["TRUE"], "Sepsis"), replicate(table(grepl("N", colnames(aggdata)))["TRUE"], "Vascular")); 
  bplot.obj$Outcome <- colNew
  bplot.obj <- melt(bplot.obj, id.var = "Outcome")
  
  Outcome_forcolor <-c("Sepsis","Vascular")
  color.codes<-as.character(c("#FF0000", "#6495ED"))
  color.names<-c("red","blue")
  df2=data.frame(Outcome_forcolor, color.codes, color.names)
  colnames(df2)[1] <- "Outcome"
  
  df <-merge(bplot.obj,df2, by=("Outcome"), all.x=TRUE, all.y=TRUE)
  
  setEPS() 
  postscript("Boxplot_immune.eps")
  ggplot(data = df, aes(x = variable, y = value, fill = Outcome)) + geom_boxplot(size = 0.75, color = "black", outlier.color = "black", outlier.shape = 16, outlier.size = 1, outlier.stroke = 2) + 
    scale_fill_manual(breaks = df$Outcome, values = unique(as.character(df$color.codes))) + scale_x_discrete(name = "Immune cell type") + scale_y_continuous(name = "Average expression of cell transcript", breaks = seq(4, 7.5, 0.25), limits = c(4, 7.5)) +
    theme_bw() + theme(axis.text = element_text(size = 11), axis.title = element_text(size=12, face = "bold")) + theme(legend.position="top", legend.title=element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  dev.off()
}


dir_reg_plt <- function(total){
  brplot <- total[,c("Fold.Change", "Cell.Specificity")]
  temp <- data.frame(Cell.Specificity = brplot$Cell.Specificity, down.regulated = (brplot$Fold.Change<0), up.regulated = (brplot$Fold.Change>0))
  temp_1 <- aggregate(down.regulated~Cell.Specificity, temp, sum)
  temp_2 <- aggregate(up.regulated~Cell.Specificity, temp, sum)
  final_binary <- merge(temp_1,temp_2,by="Cell.Specificity")
  
  setEPS()
  postscript("Barplot_immune.eps")
  ggplot(final_binary, aes(Cell.Specificity)) + 
    geom_bar(width = 0.2)+geom_bar(aes(y=-down.regulated), stat = "identity", fill = "cornflowerblue", width = 0.3) + 
    geom_bar(aes(y=up.regulated), stat = "identity", fill = "indianred", width = 0.3) + geom_hline(yintercept = 0, size = 2, color = "grey") +  #+ scale_y_continuous(labels = commapos)
    scale_x_discrete(name = "Immune cell type") + scale_y_continuous(name = "Number of probes in transcript", breaks = seq(-65, 65, 5), limits = c(-65, 65)) +
    theme_bw() + theme(axis.text = element_text(size = 11), axis.title = element_text(size=12, face = "bold")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  dev.off()
}