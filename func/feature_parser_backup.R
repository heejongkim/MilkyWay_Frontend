# 2017/2/21 - initial version: leave only unique peptides
# Hee Jong Kim
feature_parser <- function(protein_name, fasta_df, psm_table, exp_design_for_pViz){
    #protein_name <- "sp|Q9FKS0|UKL1_ARATH" #as.character(protein_saint_table()[input$protein_saint_table_row_last_clicked,]$Protein)
    protein_name <- escapeRegex(protein_name)
    protein_sequence <- as.character(fasta_df[grepl(protein_name, fasta_df$id),]$Sequence)
    filteredByProtein_psm_table <- psm_table[grepl(protein_name, psm_table$`protein id`), c("file", "file_idx", "scan", "charge", "spectrum precursor m/z", "percolator q-value", "protein id", "sequence", "proteinacc_start_stop_pre_post_;")]
    # cSplit? splited_each_chrom <- cSplit(filtered_chrom, c("Times","Intensities"), ",", "long")
    filteredByProtein_psm_table <- cSplit(filteredByProtein_psm_table, c("protein id","proteinacc_start_stop_pre_post_;"), c(",",";"), "long")
    filteredByProtein_psm_table <- filteredByProtein_psm_table[grepl(protein_name, filteredByProtein_psm_table$`protein id`), ]
    filteredByProtein_psm_table$`proteinacc_start_stop_pre_post_;` <- gsub(paste0(protein_name,'_'), '', filteredByProtein_psm_table$`proteinacc_start_stop_pre_post_;`)
    filteredByProtein_psm_table <- cSplit(filteredByProtein_psm_table, "proteinacc_start_stop_pre_post_;", "_", fixed=FALSE)
    setnames(filteredByProtein_psm_table,
             old = c("proteinacc_start_stop_pre_post_;_1","proteinacc_start_stop_pre_post_;_2","proteinacc_start_stop_pre_post_;_3","proteinacc_start_stop_pre_post_;_4"),
             new = c("x", "y", "preAA", "postAA"))
    # making duplicated columns for pViz syntax match: file to category & sequence to text
    filteredByProtein_psm_table$category <- filteredByProtein_psm_table$file
    filteredByProtein_psm_table$description <- filteredByProtein_psm_table$sequence
    # Extract only columns to pass onto pViz
    pVizFeature_df <- filteredByProtein_psm_table[,c("category", "file_idx", "description", "x", "y")]
    # set the type of annotation format
    #pVizFeature_df$type <- rep("bar", nrow(pVizFeature_df))
    # drop .mzML in file name
    pVizFeature_df$category <- gsub('.mzML', '', pVizFeature_df$category)
    # Merge pViz_feature dataframe with experimental design dataframe to fill the no id runs (potentially sort)
    pVizFeature_df <- merge(x=exp_design_for_pViz, y=pVizFeature_df, by="category", all=TRUE)
    pVizFeature_df <- data.frame(lapply(pVizFeature_df, as.character), stringsAsFactors=FALSE)
    
    # Clean up for JSON format
    pVizFeature_df$x[is.na(pVizFeature_df$x)] <- "0"
    pVizFeature_df$y[is.na(pVizFeature_df$y)] <- "0"
    #pVizFeature_df$type[is.na(pVizFeature_df$type)] <- "circle"
    pVizFeature_df$description[is.na(pVizFeature_df$description)] <- "No Identification"
    pVizFeature_df$`Biological Condition` <- NULL
    
    # Clean Up 2nd Round - By condition instead of each run
    pVizFeature_df$category <- pVizFeature_df$`Biological.Condition`
    pVizFeature_df$`Biological.Condition` <- NULL
    pVizFeature_df$file_idx <- NULL
    
    # Remove duplicates
    setkey(pVizFeature_df, NULL)
    pVizFeature_df <- unique(pVizFeature_df)
    
    return(pVizFeature_df)
  }