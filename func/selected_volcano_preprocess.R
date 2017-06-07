# Hee Jong Kim
########### volcano plot dataframe preprocessing function
selected_volcano_preprocess <- function(MSstats_comparison_df, filtered_MSstats_comparison_df, pvalue_cutoff, log2FC_cutoff){
  MSstats_comparison <- filtered_MSstats_comparison_df
  # remove completeMissing rows
  MSstats_comparison$issue <- as.character(MSstats_comparison$issue)
  MSstats_comparison[is.na(MSstats_comparison$issue),"issue"] <- ""
  MSstats_comparison <- MSstats_comparison[MSstats_comparison$issue!="completeMissing",]
  
  logFC <- colnames(MSstats_comparison)[grepl("FC", colnames(MSstats_comparison))]
  comparison_volcano_dataset <- MSstats_comparison[c("Protein", logFC, "adj.pvalue", "Label")]

  ##### INF to finite value swap part -- disabled for now
  # # Inf value handling: set biggest finite value + 1
  # logFCvector <- MSstats_comparison_df[,logFC]
  # finiteMax <- max(logFCvector[logFCvector!=max(logFCvector)]) + 1
  # finiteMin <- min(logFCvector[logFCvector!=min(logFCvector)]) - 1

  # # Replace Inf and -Inf with finiteMax/Min
  # comparison_volcano_dataset[is.infinite(comparison_volcano_dataset[,logFC]) & sign(comparison_volcano_dataset[,logFC])==1, logFC] <- finiteMax
  # comparison_volcano_dataset[is.infinite(comparison_volcano_dataset[,logFC]) & sign(comparison_volcano_dataset[,logFC])==-1, logFC] <- finiteMin

  # # adj.pvalue = 0 to the lowest adj.pvalue 
  # adjPvalueVector <- MSstats_comparison_df[, "adj.pvalue"]
  # adjPvalueVector <- adjPvalueVector[!is.na(adjPvalueVector)]
  # finiteMinAdjPvalue <- min(adjPvalueVector[adjPvalueVector!=min(adjPvalueVector)]) / 2
  # #comparison_volcano_dataset[comparison_volcano_dataset[,"adj.pvalue"]==0, "adj.pvalue"] <- finiteMinAdjPvalue
  # comparison_volcano_dataset[,"adj.pvalue"] <- ifelse(!is.na(comparison_volcano_dataset[,"adj.pvalue"]) & comparison_volcano_dataset[,"adj.pvalue"] == 0, finiteMinAdjPvalue, comparison_volcano_dataset[,"adj.pvalue"])

  comparison_volcano_dataset['group'] <- "NotSignificant"
  comparison_volcano_dataset[which(comparison_volcano_dataset['adj.pvalue']<pvalue_cutoff & abs(comparison_volcano_dataset[logFC])<log2FC_cutoff), "group"] <- "Significant"
  comparison_volcano_dataset[which(comparison_volcano_dataset['adj.pvalue']>pvalue_cutoff & abs(comparison_volcano_dataset[logFC])>log2FC_cutoff), "group"] <- "FoldChange"
  comparison_volcano_dataset[which(comparison_volcano_dataset['adj.pvalue']<pvalue_cutoff & comparison_volcano_dataset[logFC]>log2FC_cutoff), "group"] <- "SignificantANDPositiveFoldChange"
  comparison_volcano_dataset[which(comparison_volcano_dataset['adj.pvalue']<pvalue_cutoff & comparison_volcano_dataset[logFC]< -log2FC_cutoff), "group"] <- "SignificantANDNegativeFoldChange"
  return(comparison_volcano_dataset)
}
########### /volcano plot dataframe preprocessing function
