# Hee Jong Kim
##################### Functions unique for this excel exporter ############################################
########## Add basic info with the table
add_sheet_with_info <- function(wb, 
                                sheet, 
                                header, 
                                description, 
                                table){
  # Add Header
  xlsx.addHeader(wb, sheet, value = header, level = 1, color = "black", underline = 1)
  xlsx.addLineBreak(sheet, 1)
  # Add paragraph
  xlsx.addParagraph(wb, sheet,value = description, isItalic=TRUE, colSpan=5, 
                    rowSpan=4, fontColor="darkgray", fontSize=14)
  xlsx.addLineBreak(sheet, 3)
  # Add Table
  xlsx.addTable(wb, sheet, table, col.names=TRUE, row.names=TRUE, startCol = 2, rowFill=c("#CCCCCC","#FFFFFF"))
}
########## /Add basic info with the table

########### Conditional formatting cells
conditional_formatting_cells <- function(wb, 
                                         sheetname, 
                                         pvalue_cutoff, 
                                         log2FC_cutoff, 
                                         quant_table){
  # http://stackoverflow.com/questions/21618556/export-data-frames-to-excel-via-xlsx-with-conditional-formatting
  # Color Setup
  yellow_fillcolor = "#FFE892"
  yellow_fillobj <- Fill(foregroundColor=yellow_fillcolor)
  yellow_cs <- CellStyle(wb, fill=yellow_fillobj) 
  red_fillcolor = "#febfc7"
  red_fillobj <- Fill(foregroundColor=red_fillcolor)              # create fill object
  red_cs <- CellStyle(wb, fill=red_fillobj)                       # create cell style
  blue_fillcolor = "#7abad1"
  blue_fillobj <- Fill(foregroundColor=blue_fillcolor)            # create fill object
  blue_cs <- CellStyle(wb, fill=blue_fillobj)                     # create cell style
  sheets <- getSheets(wb)                                         # get all sheets
  sheet <- sheets[[sheetname]]                                    # get specific sheet
  rows <- getRows(sheet, rowIndex=12:(nrow(quant_table)+11))      # get rows
  
  # 1st row is headers
  pvalue_column_index <- list()
  log2FC_column_index <- list()
  for (i in 1:((length(quant_table)-1) / 2)){
    pvalue_column_index[i] <- 4+2*(i-1)
    log2FC_column_index[i] <- 5+2*(i-1)
  }
  pvalue_cells <- getCells(rows, colIndex = unlist(pvalue_column_index))       # get cells
  log2FC_cells <-getCells(rows, colIndex = unlist(log2FC_column_index))
  
  ##### Add red color on significant pvalue cells
  pvalue_values <- lapply(pvalue_cells, getCellValue) # extract the values
  # find cells meeting conditional criteria :: pvalue
  highlight <- "pvalue"
  for (i in names(pvalue_values)) {
    x <- as.numeric(pvalue_values[i])
    if (x <= pvalue_cutoff & !is.na(x)) {
      highlight <- c(highlight, i)
    }    
  }
  highlight <- highlight[-1]
  lapply(names(pvalue_cells[highlight]),
         function(ii)setCellStyle(pvalue_cells[[ii]],yellow_cs))
  
  ##### Add red color on high positive fold change cells
  log2FC_values <- lapply(log2FC_cells, getCellValue) # extract the values
  # find cells meeting conditional criteria :: log2FC
  highlight <- "log2FC"
  for (i in names(log2FC_values)) {
    x <- as.numeric(log2FC_values[i])
    if (x >= log2FC_cutoff & !is.na(x)) {
      highlight <- c(highlight, i)
    }
  }
  highlight <- highlight[-1]
  lapply(names(log2FC_cells[highlight]),
         function(ii)setCellStyle(log2FC_cells[[ii]],red_cs))
  
  ##### Add blue color on high negative fold change cells
  log2FC_values <- lapply(log2FC_cells, getCellValue) # extract the values
  # find cells meeting conditional criteria :: log2FC
  highlight <- "log2FC"
  for (i in names(log2FC_values)) {
    x <- as.numeric(log2FC_values[i])
    if (x <= -log2FC_cutoff & !is.na(x)) {
      highlight <- c(highlight, i)
    }
  }
  highlight <- highlight[-1]
  lapply(names(log2FC_cells[highlight]),
         function(ii)setCellStyle(log2FC_cells[[ii]],blue_cs))
}
########### /Conditional formatting cells

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


##################### /Functions unique for this excel exporter ############################################

excel_report_generator <- function(filename,
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
                                   quantification_csv
                                   ){
  # Initiate new Excel workbook
  wb <- createWorkbook(type="xlsx")

  ##################### Exp_design TAB
  # Create sheet for Experiment design table
  sheet <- createSheet(wb, sheetName = "Exp_design")
  experiment_design$`Test or Control` <- as.character(experiment_design$`Test or Control`)
  # Experimental design info and table
  add_sheet_with_info(wb,
                      sheet = sheet,
                      header = "Experimental Design for Statistical Analysis of Mass Spectrometry Dataset",
                      description = paste("Brief Description of \n",
                                          "this table can be here. Along with \n",
                                          "statistical analysis tool's publication information", sep=""),
                      table = experiment_design[order(experiment_design["Test or Control"], experiment_design["BioReplicate"], experiment_design["Original File Name"]), ])
  # Experiment Diagram Insertion
  exp_diagram(experiment_design) %>%
    export_graph(file_name = "exp_diagram.png",
                 file_type = "PNG")
  nrows<-length(getRows(sheet))   # Count the table row
  addPicture("exp_diagram.png", sheet, startColumn = 2, startRow = nrows + 2)   # Insert the diagram at the bottom
  xlsx.addLineBreak(sheet, 2)   # Add Break line
  ####################################### /Exp_design TAB

  ####################################### Quant TAB
  # Generate quant_table
  # Create pivot tables
  comparison_pvalue_table <- dcast(comparison_csv, Protein ~ Label, value.var = "adj.pvalue")
  colnames(comparison_pvalue_table) <- paste(colnames(comparison_pvalue_table), "adj.pvalue", sep="_")
  comparison_log2FC_table <- dcast(comparison_csv, Protein ~ Label, value.var = colnames(comparison_csv)[grepl("FC", colnames(comparison_csv))])
  colnames(comparison_log2FC_table) <- paste(colnames(comparison_log2FC_table), colnames(comparison_csv)[grepl("FC", colnames(comparison_csv))], sep="_")
  # Merge pivot tables
  FC_name = paste0("Protein_", colnames(comparison_csv)[grepl("FC", colnames(comparison_csv))])
  quant_table <- merge(x=comparison_pvalue_table, y=comparison_log2FC_table, by.x="Protein_adj.pvalue", by.y=FC_name, all = TRUE)
  colnames(quant_table)[colnames(quant_table) == 'Protein_adj.pvalue'] <- "Protein"
  # Reorder the columns by conditions
  setcolorder(quant_table, c("Protein",sort(colnames(quant_table)[2:length(colnames(quant_table))])))
  # Add Fasta Description
  quant_table <- merge(x=quant_table, y=fasta_df, by.x="Protein", by.y="id", all.x=TRUE)
  quant_table$Sequence <- NULL
  # Create sheet
  sheet <- createSheet(wb, sheetName = "Quant_table")
  # pvalue & log2FC pivot table from MSstats comparison_csv
  add_sheet_with_info(wb,
                      sheet = sheet,
                      header = "Pivot Tables from MSStats",
                      description = paste("Brief Description of \n",
                                          "this merged pivot table can be here. Along with \n",
                                          "statistical analysis tool's publication information", sep=""),
                      table = quant_table)
  # Add filters
  addAutoFilter(sheet, paste(paste0(LETTERS[3],11),paste0(LETTERS[ncol(quant_table)+2],nrow(quant_table)+11), sep=":"))
  # conditional formatting to highlight significant cells
  conditional_formatting_cells(wb,
                               sheetname = "Quant_table",
                               pvalue_cutoff = pvalue_cutoff,
                               log2FC_cutoff = log2FC_cutoff,
                               quant_table = quant_table)
  ####################################### /Quant TAB

  ####################################### ID TAB
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
  # Create sheet
  sheet <- createSheet(wb, sheetName = "ID_table")
  # Saint table from saint_table
  add_sheet_with_info(wb,
                      sheet = sheet,
                      header = "SAINTexpress Analysis and SpC Table",
                      description = paste("Brief Description of \n",
                                          "this table can be here. Along with \n",
                                          "statistical analysis tool's publication information", sep=""),
                      table = id_table)
  # Add filters
  addAutoFilter(sheet, paste(paste0(LETTERS[3],11),paste0(LETTERS[ncol(id_table)+2],nrow(id_table)+11), sep=":"))
  ####################################### /ID TAB

  ####################################### Figure Tab
  # Create sheet
  sheet <- createSheet(wb, sheetName = "Figures")
  # Create multiple volcano plots
  volcano_file_list <- multiple_volcano_plots(comparison_csv = comparison_csv,
                                              pvalue_cutoff = pvalue_cutoff,
                                              log2FC_cutoff = log2FC_cutoff)
  # Add volcano plots into the sheet
  for(i in 1:length(volcano_file_list)){
    addPicture(volcano_file_list[[i]], sheet, scale = 0.3, startColumn = 4 + 10 * (i-1), startRow = 5)
  }

  # Create UpSet Plots
  psm_upsetR(experiment_design,
             psm_table,
             peptide_q_value = peptide_q_value,
             protein_q_value = protein_q_value,
             peptidoform_img_filename = "peptidoform_upset.png",
             protein_img_filename = "protein_upset.png")
  # Add peptidoform upset & protein upset plots
  addPicture("peptidoform_upset.png", sheet, scale = 0.3, startColumn = 4 + 10 * 0, startRow = 52)
  addPicture("protein_upset.png", sheet, scale = 0.3, startColumn =  4 + 10 * 1, startRow = 52)

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
  addPicture("heatmap.png", sheet, scale = 0.3, startColumn = 4 + 10 * 0, startRow = 99)
  # Clean up
  file.remove("exp_diagram.png")
  file.remove("heatmap.png")
  file.remove(list.files(pattern = "volcano_plot_*"))
  file.remove(list.files(pattern = "*_upset.png"))
  ####################################### Figure Tab
  # Save
  saveWorkbook(wb, filename)
  return(filename)
}