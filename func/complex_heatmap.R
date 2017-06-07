########### complex heatmap plotting function
# Hee Jong Kim & William Barshop
complex_heatmap <- function(unidirectional, quantCentering, minPvalueFilter, experimental_design_file, comparison_data, quantification_data,manual_k){
  #Grab list of conditions used as controls...
  controls<-unique(experimental_design_file[experimental_design_file$`Test or Control`=="C","Biological Condition"])
  controls<-gsub("\\+","\\\\\\+",controls)
  experimental_design_file$Sample<-paste(experimental_design_file[,"Biological Condition"],experimental_design_file[,"BioReplicate"],sep="_")
  # Set protein name column as row name
  rownames(quantification_data) <- quantification_data$Protein
  quantification_data$Protein <- NULL
  quantification_data$X <- NULL
  ##########################################################################
  ### Create a pivot dataframe from log2FC with at least one condition is positive
  ##########################################################################
  # make pivot table from comparison_data Protein row and Label column with value from log2FC column
  pivot_log2FC <- dcast(comparison_data, Protein ~ Label, value.var = colnames(comparison_data)[grepl("FC", colnames(comparison_data))])
  # set protein name as row name
  rownames(pivot_log2FC) <- pivot_log2FC$Protein
  # remove non row name - protein name column to avoid duplication
  pivot_log2FC$Protein <- NULL
  # make new column with the maximum value in each row
  pivot_log2FC$max <- apply(pivot_log2FC,1,max)
  # extract the subset of dataframe by log2FC > 0.0 subset_pivot_log2FC <- pivot_log2FC[pivot_log2FC$max > 0.0, ]
  subset_pivot_log2FC <- pivot_log2FC
  # remove "max" column
  subset_pivot_log2FC$max <- NULL
  # replace log2FC value to binary -1/1 values ---- at this point, negative direction is set to zero if unidirectional
  if(unidirectional){
    subset_pivot_log2FC[subset_pivot_log2FC<0] <- 0
  } else {
    subset_pivot_log2FC[subset_pivot_log2FC<0] <- -1
  }
  subset_pivot_log2FC[subset_pivot_log2FC>0] <- 1
  ##########################################################################
  ### Create a pivot dataframe from agj. p-value with at least one condition is significant (<=0.05)
  ##########################################################################
  # make pivot table from comparison_data Protein row and Label column with value from adj.pvalue column
  pivot_adj_pvalue <- dcast(comparison_data, Protein ~ Label, value.var = "adj.pvalue")
  # set protein name as row name
  rownames(pivot_adj_pvalue) <- pivot_adj_pvalue$Protein
  # remove non row name - protein name column to avoid duplication
  pivot_adj_pvalue$Protein <- NULL
  # make new column with the minimum value in each row
  pivot_adj_pvalue$min <- apply(pivot_adj_pvalue,1,min)
  # extract the subset of dataframe by adjusted p-value <= 0.05
  subset_pivot_adj_pvalue <- pivot_adj_pvalue[pivot_adj_pvalue$min <= minPvalueFilter, ]
  # remove "min" column
  subset_pivot_adj_pvalue$min <- NULL
  ##########################################################################
  ##########################################################################
  ### Create directional probability table
  ##########################################################################
  # obtain column names
  exp_column_name <- colnames(subset_pivot_adj_pvalue)
  # inner join "subset_pivot_adj_pvalue" and "subset_pivot_log2FC"
  joined_pivot_adj_pvalue_log2FC <- merge(x=subset_pivot_adj_pvalue, y=subset_pivot_log2FC, by="row.names")
  # Set the protein name column as row name again!
  rownames(joined_pivot_adj_pvalue_log2FC) <- joined_pivot_adj_pvalue_log2FC$Row.names
  joined_pivot_adj_pvalue_log2FC$Row.names <- NULL
  #Adjust to directional probabilities or unidirectional "Enrichment Probabilities"
  for (each_column_name in exp_column_name){
    adj_pvalue_col_name <- paste(each_column_name, ".x", sep="")
    log2FC_direction_col_name <- paste(each_column_name, ".y", sep="")
    joined_pivot_adj_pvalue_log2FC[,each_column_name] <- (1-joined_pivot_adj_pvalue_log2FC[,adj_pvalue_col_name]) * joined_pivot_adj_pvalue_log2FC[,log2FC_direction_col_name]
  }
  # Drop all .x and .y columns
  comparison_directional_p <- joined_pivot_adj_pvalue_log2FC[, names(joined_pivot_adj_pvalue_log2FC) %in% exp_column_name, drop = F]
  ## We're only interested in those enriched in at least 1 sample, if unidirectional, otherwise, this filtering is already taken care of...
  if(unidirectional){
    comparison_directional_p$max <- apply(comparison_directional_p,1,max)
    comparison_directional_p <- comparison_directional_p[comparison_directional_p$max >= (1-minPvalueFilter), ]
    comparison_directional_p$max <- NULL
  }
  ##########################################################################
  ### Merge to distill the rows down for heatmap
  ##########################################################################
  # by join the qaunt_data with adj_pvalue to match the rows (filtering)
  adj_pvalue_quant_joint <- merge(x = comparison_directional_p, y = quantification_data, by="row.names", all.x=TRUE)
  # drop the temporary adj_pvalue columns
  subset_quant <- adj_pvalue_quant_joint[, ! names(adj_pvalue_quant_joint) %in% colnames(comparison_directional_p), drop = F]
  # Set the protein name column as row name again!
  rownames(subset_quant) <- subset_quant$Row.names
  subset_quant$Row.names <- NULL
  ##########################################################################
  ###### RE CORRECT THE DASHES WHICH WERE REMOVED...
  names(subset_quant) <- gsub(x = names(subset_quant),pattern = "\\.",
                              replacement = "-")
  #Cleanup Merge Column junk...
  subset_quant$Var.2<-NULL
  subset_quant$`Var.2`<-NULL
  subset_quant$`Var-2`<-NULL
  #print(subset_quant)
  
  ##### Let's work on setting colors correctly...
  quantnorm<-subset_quant
  #print(quantnorm)
  if(quantCentering=="CONTROL"){
    #print(paste(controls,sep='|')) print(names(quantnorm)) print(grep(paste(controls,sep='|'),names(quantnorm),value = TRUE)) print(names(quantnorm)) print(quantnorm[,grep(paste(controls,sep='|'),value = TRUE)])
    if(length(grep(paste(controls,sep='|'),names(quantnorm),value = TRUE))>1){
      mean_std<-rowMeans(quantnorm[,grep(paste(controls,sep='|'),names(quantnorm),value = TRUE)],na.rm = TRUE)
    }else{
      mean_std<-quantnorm[,grep(paste(controls,sep='|'),names(quantnorm),value = TRUE)]
    }
    quantnorm<-quantnorm-mean_std
    intensityGraphName<-"log2(intensity) minus\nmean control value"
  } else if(quantCentering=="ALL"){
    mean_all<-rowMeans(quantnorm)
    quantnorm<-quantnorm-mean_all
    intensityGraphName<-"log2(intensity)\nminus mean"
  } else if(quantCentering=="NONE"){
    print("not centering intensity values")
    intensityGraphName<-"log2(intensity)"
  }
  #If unidirectional, we're working from 0-->1, otherwise -1-->0-->1
  if(unidirectional){
    probabilityGraphName<-"Enrichment Probability"
    probability_colorRamp<-colorRamp2(c(0,1),c("#000000","#FF0000"))
  } else {
    probabilityGraphName<-"Directional Probability"
    probability_colorRamp<-colorRamp2(c(-1,0,1),c("#0000FF","#000000","#FF0000"))
  }
  #Quant colors will range from blue at minimum, to black at 0, to red at maximum print(paste("min",min(as.numeric(quantnorm)),"max",max(as.numeric(quantnorm)))) print(quantnorm) print(min(as.numeric(quantnorm))) 
  #print(max(as.numeric(quantnorm)))
  quant_colorRamp<-colorRamp2(c(min(quantnorm,na.rm=TRUE),0,max(quantnorm,na.rm=TRUE)),c("#0000FF","#000000","#FF0000"))
  ####Below here is the color mapping information...
  color_mapping_df<-data.frame(Condition=unique(experimental_design_file$`Biological Condition`))
  colors_to_use<-rainbow(length(unique(color_mapping_df$Condition)))
  for(i in seq(1,length(unique(color_mapping_df$Condition)))){
    color_mapping_df[i,"color"]<-colors_to_use[i]
  }
  mapping_df<-data.frame(Condition=experimental_design_file$`Biological Condition`,Run=experimental_design_file$Sample)
  groups<-mapping_df[match(colnames(quantnorm),mapping_df$Run),1,drop=F]
  strip <- function (x) gsub("^\\s+|\\s+$", "", x)
  colnames(comparison_directional_p)<-sapply(strsplit(colnames(comparison_directional_p),"vs.",1),"[",1)
  names(comparison_directional_p) <- strip(names(comparison_directional_p))
  directional_prob_colors<-mapvalues(colnames(comparison_directional_p),from=color_mapping_df$Condition,to=color_mapping_df$color)
  quant_colors<-mapvalues(as.vector(groups['Condition'])[[1]],from=as.vector(color_mapping_df$Condition),to=as.vector(color_mapping_df$color))
  prob_color_list<-character()
  for(j in seq(1,length(colnames(comparison_directional_p)))){
    prob_color_list<-c(prob_color_list,toString(directional_prob_colors[j]))
    names(prob_color_list)[length(prob_color_list)]<-toString(colnames(comparison_directional_p)[j])
  }
  quant_color_list<-character()
  for(j in seq(1,length(colnames(quantnorm)))){
    quant_color_list<-c(quant_color_list,toString(quant_colors[j]))
    names(quant_color_list)[length(quant_color_list)]<-toString(groups[j,'Condition'])
  }
  #Now we'll add back the controls which might be missing from the probability color mappings!
  for(j in seq(1,length(quant_color_list[controls]))){
    if(!quant_color_list[j] %in% prob_color_list){
      prob_color_list<-append(prob_color_list,quant_color_list[j])
      names(prob_color_list)[length(prob_color_list)]<-toString(names(quant_color_list[j]))
    }
  }
  ha = HeatmapAnnotation(df=groups['Condition'],col=list(Condition=quant_color_list),show_legend = FALSE)
  ha2 = HeatmapAnnotation(df=data.frame(Condition=colnames(comparison_directional_p)),col=list(Condition=prob_color_list))
  # Setup heatmaps
  max_num_clusters<-15 #This is the default
  
  
  
  
  
  
  #### If we take in the input to handle cluster size manually... do it manually.
  if(manual_k<0){
    
    if (length(comparison_directional_p[,1])<15 & length(comparison_directional_p[,1])>4){
      max_num_clusters<-4
      nbc<-NbClust(comparison_directional_p,method="kmeans",max.nc=max_num_clusters)
    } else if (length(comparison_directional_p[,1])<=4){
      nbc<-1
    } else {
      nbc<-NbClust(comparison_directional_p,method="kmeans",max.nc=max_num_clusters)
    }
    optimal_k<-as.integer(names(sort(table(nbc$Best.nc[1,]),decreasing=TRUE)[1])[1])
    if (is.null(optimal_k)){
      optimal_k<-1
    }
  } else {
    optimal_k<-manual_k
  }
  heatmap_adj_pvalue <- Heatmap(comparison_directional_p,
                                km=optimal_k,
                                row_title_side="left",
                                row_names_side = "right",
                                row_dend_side="left",
                                row_dend_width = unit(20, "mm"),
                                column_title = probabilityGraphName,
                                name = probabilityGraphName,
                                gap = unit(5, "mm"),
                                col=probability_colorRamp,
                                top_annotation = ha2)
  heatmap_quant <- Heatmap(quantnorm,
                           column_title = intensityGraphName,
                           name = "Intensity",
                           cluster_columns = FALSE,
                           show_row_names = FALSE,
                           col=quant_colorRamp,
                           top_annotation = ha)
  # Generate multiple heatmaps
  return(heatmap_adj_pvalue + heatmap_quant)
}
########### /complex heatmap plotting function