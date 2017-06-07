########### Generate Experiment table, Quant table, and Identificaiton table
# Hee Jong Kim
generate_tables_qual <- function(experiment_design, 
                            saint_table, 
                            spc_table,
                            fasta_df){
  ##### Exp_design TAB
  experiment_design$`Test or Control` <- as.character(experiment_design$`Test or Control`)
  exp_design_table <- experiment_design[order(experiment_design["Test or Control"], experiment_design["BioReplicate"], experiment_design["Original File Name"]), ]
    
  ##### ID TAB
  # Add suffix on colnames
  # Check if SAINT table exists or not
  if(!any(grepl("FIDO", colnames(saint_table)))){
    saint_table_excel <- saint_table
    colnames(saint_table_excel) <- paste(colnames(saint_table_excel), "SAINT", sep="_")
    colnames(saint_table_excel)[1] <- "Prey"
  }
  spc_table_excel <- spc_table
  colnames(spc_table_excel) <- paste(colnames(spc_table_excel), "SpC", sep="_")
  colnames(spc_table_excel)[1] <- "ProteinID"
  # Merge saint and Spc
  if(!any(grepl("FIDO", colnames(saint_table)))){
    id_table <- merge(x=saint_table_excel, y=spc_table_excel, by.x="Prey", by.y="ProteinID", all=TRUE)
    colnames(id_table)[colnames(id_table)=="Prey"] <- "Protein"
  }else{
    id_table <- spc_table_excel
    colnames(id_table)[colnames(id_table)=="ProteinID"] <- "Protein"
  }
  
  id_table <- merge(x=id_table, y=fasta_df, by.x="Protein", by.y="id", all.x=TRUE)
  if(!any(grepl("FIDO", colnames(saint_table)))){
    id_table$description_SAINT <- NULL
  }
  id_table$Sequence <- NULL
  
  # Combine all three tables in a list
  excel_sheets <- list("Exp_Design" = exp_design_table,
                       "Ident_Table" = id_table)
  return(excel_sheets)
}


########### Create UpSet Plots from Peptidoform count and protein count per condition
# Identical function
########### /Create UpSet Plots from Peptidoform count and protein count per condition


########### Excel report generator main function
excel_report_generator_qual <- function(filename,
                                      # UpSet Plot
                                      peptide_q_value,
                                      protein_q_value,
                                      # data from RData
                                      experiment_design,
                                      fasta_df,
                                      saint_table = NA,
                                      spc_table,
                                      psm_table){

  ##### ADD all three tables
  excel_sheets <- generate_tables_qual(experiment_design, saint_table, spc_table, fasta_df)
  wb <- write.xlsx(excel_sheets, file = filename, asTable = TRUE, tableStyle = "TableStyleLight16")
  # Add style to each sheet
  setColWidths(wb, sheet = 1, cols=1:7,widths=40) # Exp_Design
  setColWidths(wb, sheet = 2, cols=1:ncol(excel_sheets$Ident_Table),widths=20) # ID Table
  ##### /ADD all three tables
  
  # Initial figures sheet setup
  addWorksheet(wb, "Figures", gridLines = FALSE)
  ##### Experiment Diagram Insertion
  exp_diagram(experiment_design) %>%
    export_graph(file_name = "exp_diagram.png",
                 file_type = "PNG")
  insertImage(wb, "Figures", "exp_diagram.png", startCol = 4, startRow = 2, width = 4000, height = 2000, units="px")
  ##### /Experiment Diagram Insertion
  
  ##### ADD UpSet Plots
  # Create UpSet Plots
  if(length(unique(psm_table$file))>1){
      psm_upsetR(experiment_design,
                 psm_table,
                 peptide_q_value = peptide_q_value,
                 protein_q_value = protein_q_value,
                 peptidoform_img_filename = "peptidoform_upset.png",
                 protein_img_filename = "protein_upset.png")
      # Add peptidoform upset & protein upset plots
      insertImage(wb, "Figures", "peptidoform_upset.png", startCol = 4 + 10 * 0, startRow = 90, width = 2100, height = 3000, units="px")
      insertImage(wb, "Figures", "protein_upset.png", startCol =  4 + 10 * 1, startRow = 90, width = 2100, height = 3000, units="px")
    }
  ##### /ADD UpSet Plots
  
  ##### FINALLY SAVE THEM ALL
  saveWorkbook(wb, filename, overwrite=TRUE)
  
  ##### Clean up all figure images
  file.remove("exp_diagram.png")
  file.remove(list.files(pattern = "*_upset.png"))
  return(filename)
}