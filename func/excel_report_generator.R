########### Generate Experiment table, Quant table, and Identificaiton table
# Hee Jong Kim
generate_tables <- function(experiment_design, 
                            comparison_csv, 
                            saint_table, 
                            spc_table,
                            outlier_table,
                            fasta_df){
  ##### Exp_design TAB
  experiment_design$`Test or Control` <- as.character(experiment_design$`Test or Control`)
  exp_design_table <- experiment_design[order(experiment_design["Test or Control"], experiment_design["BioReplicate"], experiment_design["Original File Name"]), ]
  
  ##### Quant TAB
  # Generate quant_table
  # Create pivot tables
  logFC_colname <- colnames(comparison_csv)[grepl("FC", colnames(comparison_csv))]
  comparison_pvalue_table <- dcast(comparison_csv, Protein ~ Label, value.var = "adj.pvalue")
  colnames(comparison_pvalue_table) <- paste(colnames(comparison_pvalue_table), "adj.pvalue", sep="_")
  comparison_log2FC_table <- dcast(comparison_csv, Protein ~ Label, value.var = logFC_colname)
  colnames(comparison_log2FC_table) <- paste(colnames(comparison_log2FC_table), logFC_colname, sep="_")
  # Merge pivot tables
  FC_name = paste0("Protein_", logFC_colname)
  quant_table <- merge(x=comparison_pvalue_table, y=comparison_log2FC_table, by.x="Protein_adj.pvalue", by.y=FC_name, all = TRUE)
  colnames(quant_table)[colnames(quant_table) == 'Protein_adj.pvalue'] <- "Protein"
  # Reorder the columns by conditions
  setcolorder(quant_table, c("Protein",sort(colnames(quant_table)[2:length(colnames(quant_table))])))
  # Add Fasta Description
  quant_table <- merge(x=quant_table, y=fasta_df, by.x="Protein", by.y="id", all.x=TRUE)
  quant_table$Sequence <- NULL
  
  ##### ID TAB
  # Add suffix on colnames
  saint_table_excel <- saint_table
  colnames(saint_table_excel) <- paste(colnames(saint_table_excel), "SAINT", sep="_")
  colnames(saint_table_excel)[1] <- "Prey"
  spc_table_excel <- spc_table
  colnames(spc_table_excel) <- paste(colnames(spc_table_excel), "SpC", sep="_")
  colnames(spc_table_excel)[1] <- "ProteinID"
  # Merge saint and Spc
  id_table <- merge(x=saint_table_excel, y=spc_table_excel, by.x="Prey", by.y="ProteinID", all=TRUE)
  colnames(id_table)[colnames(id_table)=="Prey"] <- "Protein"
  id_table <- merge(x=id_table, y=fasta_df, by.x="Protein", by.y="id", all.x=TRUE)
  id_table$Sequence <- NULL
  id_table$description_SAINT <- NULL
  
  # Combine all three tables in a list
  excel_sheets <- list("Exp_Design" = exp_design_table,
                       "Quant_Table" = quant_table,
                       "Ident_Table" = id_table,
                       "Outlier_Table" = outlier_table)
  return(excel_sheets)
}

########### Create multiple volcano plots
multiple_volcano_plots <- function(comparison_csv, 
                                   pvalue_cutoff, 
                                   log2FC_cutoff){
  # Make list of variable names to loop over.
  label_list = unique(comparison_csv$Label)
  
  # Make plots.
  file_list = list()
  for (i in 1:length(unique(comparison_csv$Label))) {
    volcanoData <- volcano_preprocess(comparison_csv, as.character(label_list[[i]]), pvalue_cutoff, log2FC_cutoff)
    # Handling log2FC vs logFC column name here due to the quirk on ggplot colname notation
    if ("logFC" %in% colnames(volcanoData)){
      volcanoData$log2FC <- volcanoData$logFC
    }
    # manual coloring group
    volcano_group.colors <- c(NotSignificant = "#999999", 
                              Significant = "#d6d6d6", 
                              FoldChange = "#7f7f7f", 
                              SignificantANDNegativeFoldChange = "#0072B2", 
                              SignificantANDPositiveFoldChange = "#D55E00")
    
    p <- ggplot(volcanoData,
                aes_string(x=colnames(volcanoData)[grepl("FC", colnames(volcanoData))], 
                           y="-log2(adj.pvalue)",
                           key="Protein",
                           color="group")) + 
      geom_point() +
      scale_color_manual(values=volcano_group.colors) +
      #geom_text_repel(data=subset(volcanoData, abs(log2FC)>=log2FC_cutoff & adj.pvalue<=pvalue_cutoff), 
      #                aes(label=Protein)) +
      ggtitle(as.character(label_list[[i]])) +
      theme(plot.title = element_text(hjust=0.5, lineheight=.8, face="bold", size=12),
            legend.position="bottom", legend.direction="vertical")
    file_name = paste("volcano_plot_", i, ".png", sep="")
    ggsave(file_name, width = 7, height = 10, dpi = 300)
    file_list[[i]] = file_name
  }
  
  return(file_list)
}
########### /Create multiple volcano plots
########### Create UpSet Plots from Peptidoform count and protein count per condition
psm_upsetR <- function(experiment_design, 
                       psm_table, 
                       peptide_q_value, 
                       protein_q_value, 
                       peptidoform_img_filename, 
                       protein_img_filename){
  # Global preprocessing by adding .mzML and count the number of conditions
  experiment_design$file <- paste0(experiment_design$`Original File Name`, ".mzML")
  number_of_sets = as.numeric(length(unique(experiment_design$`Biological Condition`)))
  
  # peptidoform per condition UpsetR
  peptidoform_for_upsetR <- merge(x = experiment_design[,c("file", "Biological Condition")],
                                  y = psm_table[psm_table$`percolator q-value`<=peptide_q_value, c("file", "sequence")],
                                  by = "file",
                                  all = TRUE)
  peptidoform_for_upsetR$file <- NULL
  condition_peptidoforms <- dcast(peptidoform_for_upsetR, sequence ~ `Biological Condition`)
  condition_peptidoforms$sequence <- NULL
  condition_peptidoforms[condition_peptidoforms > 0] <- 1
  
  png(peptidoform_img_filename, width = 7, height = 10, units = "in", res = 300)
  upset(condition_peptidoforms, 
        nsets = number_of_sets,
        mainbar.y.label = "Peptide Count",
        sets.x.label = "Peptide per Condition", 
        sets.bar.color = "#56B4E9",
        nintersects  = NA,
        order.by = "freq", 
        empty.intersections = "on",
        keep.order = TRUE)
  dev.off()
  
  # protein per condition UpsetR
  psm_table$`protein q-values` <- as.numeric(as.character(psm_table$`protein q-values`)) # factor to numeric
  proteins_for_upsetR <- merge(x = experiment_design[,c("file", "Biological Condition")],
                               y = psm_table[psm_table$`protein q-values` <= protein_q_value, c("file", "protein id")],
                               by = "file",
                               all = TRUE)
  proteins_for_upsetR$file <- NULL
  condition_proteins <- dcast(proteins_for_upsetR, `protein id` ~ `Biological Condition`)
  condition_proteins$`protein id` <- NULL
  condition_proteins[condition_proteins > 0] <- 1
  
  png(protein_img_filename, width = 7, height = 10, units = "in", res = 300)
  upset(condition_proteins, 
        nsets = number_of_sets,
        mainbar.y.label = "Protein Count",
        sets.x.label = "Protein per Condition", 
        sets.bar.color = "#56B4E9",
        nintersects  = NA,
        order.by = "freq", 
        empty.intersections = "on",
        keep.order = TRUE)
  dev.off()
}
########### /Create UpSet Plots from Peptidoform count and protein count per condition


########### Excel report generator main function
excel_report_generator_v2 <- function(filename,
                                      # Volcano Plot, Conditional Format for Quant Table
                                      pvalue_cutoff,
                                      log2FC_cutoff,
                                      # UpSet Plot
                                      peptide_q_value,
                                      protein_q_value,
                                      # For heatmap
                                      unidirectional,
                                      quantCentering,
                                      minPvalueFilter,
                                      cluster_k,
                                      # data from RData
                                      experiment_design,
                                      comparison_csv,
                                      fasta_df,
                                      saint_table,
                                      spc_table,
                                      psm_table,
                                      quantification_csv,
                                      outlier_table){

  # Extract comparison_csv without outliers
  logFC_colname <- colnames(comparison_csv)[grepl("FC", colnames(comparison_csv))]
  comparison_csv_no_outliers <- comparison_csv[!is.infinite(comparison_csv[,logFC_colname]) & !is.na(comparison_csv[,"adj.pvalue"]),]     # Remove NA or Inf/-Inf rows for Quant tab
  ##### ADD all three tables
  excel_sheets <- generate_tables(experiment_design, comparison_csv_no_outliers, saint_table, spc_table, outlier_table, fasta_df)
  wb <- write.xlsx(excel_sheets, file = filename, asTable = TRUE, tableStyle = "TableStyleLight16")
  # Add style to each sheet
  setColWidths(wb, sheet = 1, cols=1:7,widths=40) # Exp_Design
  setColWidths(wb, sheet = 2, cols=1:ncol(excel_sheets$Quant_Table),widths=20) # Quant Table
  setColWidths(wb, sheet = 3, cols=1:ncol(excel_sheets$Ident_Table),widths=20) # ID Table
  setColWidths(wb, sheet = 4, cols=1:ncol(excel_sheets$Outlier_Table),widths=20) # Outlier Table
  ##### /ADD all three tables
  
  ##### Conditional format for Quant table sheet
  redStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
  greenStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
  yellowStyle <- createStyle(fontColour = "#9C6500", bgFill = "#FFEB9C")
  greyStyle <- createStyle(fontColour = "#2A2A2A", bgFill = "#D3D3D3")
  for (i in grep("pvalue", colnames(excel_sheets$Quant_Table))){
    conditionalFormatting(wb, sheet = 2, cols=i, rows=1:nrow(excel_sheets$Quant_Table)+1, rule=c(0,pvalue_cutoff), style = yellowStyle, type = "between")
    conditionalFormatting(wb, sheet = 2, cols=i, rows=1:nrow(excel_sheets$Quant_Table)+1, rule='=""', style = greyStyle)    # Grey for empty cell
  }
  for (i in grep("FC", colnames(excel_sheets$Quant_Table))){
    conditionalFormatting(wb, sheet = 2, cols=i, rows=1:nrow(excel_sheets$Quant_Table)+1, rule=paste0(">",log2FC_cutoff), style = redStyle, type = "expression")
    conditionalFormatting(wb, sheet = 2, cols=i, rows=1:nrow(excel_sheets$Quant_Table)+1, rule=paste0("<","-",log2FC_cutoff), style = greenStyle, type = "expression")
    conditionalFormatting(wb, sheet = 2, cols=i, rows=1:nrow(excel_sheets$Quant_Table)+1, rule='=""', style = greyStyle)    # Grey for empty cell
  }
  #####
  
  # Initial figures sheet setup
  addWorksheet(wb, "Figures", gridLines = FALSE)
  ##### Experiment Diagram Insertion
  exp_diagram(experiment_design) %>%
    export_graph(file_name = "exp_diagram.png",
                 file_type = "PNG")
  insertImage(wb, "Figures", "exp_diagram.png", startCol = 4, startRow = 2, width = 4000, height = 2000, units="px")
  ##### /Experiment Diagram Insertion
  
  ##### ADD Volcano plots
  # Create multiple volcano plots
  volcano_file_list <- multiple_volcano_plots(comparison_csv = comparison_csv,
                                              pvalue_cutoff = pvalue_cutoff,
                                              log2FC_cutoff = log2FC_cutoff)
  # Add volcano plots into the sheet
  for(i in 1:length(volcano_file_list)){
    insertImage(wb, "Figures", volcano_file_list[[i]], startCol = 4 + 10 * (i-1), startRow = 35, width = 2100, height = 3000, units="px")
  }
  ##### /ADD Volcano plots
  
  ##### ADD UpSet Plots
  # Create UpSet Plots
  psm_upsetR(experiment_design,
             psm_table,
             peptide_q_value = peptide_q_value,
             protein_q_value = protein_q_value,
             peptidoform_img_filename = "peptidoform_upset.png",
             protein_img_filename = "protein_upset.png")
  # Add peptidoform upset & protein upset plots
  insertImage(wb, "Figures", "peptidoform_upset.png", startCol = 4 + 10 * 0, startRow = 90, width = 2100, height = 3000, units="px")
  insertImage(wb, "Figures", "protein_upset.png", startCol =  4 + 10 * 1, startRow = 90, width = 2100, height = 3000, units="px")
  ##### /ADD UpSet Plots
  
  ##### ADD Heatmap
  # Generate heatmap
  png("heatmap.png", width = 7, height = 10, units = "in", res = 300)
  print(complex_heatmap(unidirectional = unidirectional,
                        quantCentering = quantCentering,
                        minPvalueFilter = minPvalueFilter,
                        experimental_design_file = experiment_design,
                        comparison_data = comparison_csv,
                        quantification_data = quantification_csv,
                        manual_k = cluster_k))
  dev.off()
  # Add heatmap
  insertImage(wb, "Figures", "heatmap.png",  startCol = 4 + 10 * 0, startRow = 140, width = 2100, height = 3000, units="px")
  ##### /ADD Heatmap
  
  ##### FINALLY SAVE THEM ALL
  saveWorkbook(wb, filename, overwrite=TRUE)
  
  ##### Clean up all figure images
  file.remove("exp_diagram.png")
  file.remove("heatmap.png")
  file.remove(list.files(pattern = "volcano_plot_*"))
  file.remove(list.files(pattern = "*_upset.png"))
  return(filename)
}