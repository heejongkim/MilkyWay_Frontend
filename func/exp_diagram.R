# Hee Jong Kim
# DiagrammeR
exp_diagram <- function(experiment_design_df){
  # Preprocessing
  exp_design <- arrange.vars(experiment_design_df, c("Test or Control"=1,
                                                     "Biological Condition"=2,
                                                     "BioReplicate"=3,
                                                     "Fractionation Group Name"=4,
                                                     "Fractionation Group ID String"=5,
                                                     "Original File Name"=6))
  exp_design <- exp_design[order(exp_design["Test or Control"], exp_design["BioReplicate"], exp_design["Original File Name"]), ]
  
  # avector <- as.character(experiment_design[['Test or Control']])
  node_vector <- c()
  node_vector <- append(node_vector, as.character(unique(exp_design[order(exp_design["Test or Control"]),"Test or Control"])))
  node_vector <- append(node_vector, as.character(unique(exp_design[,"Biological Condition"])))
  node_vector <- append(node_vector, as.character(paste("BioRep",unique(exp_design[,"BioReplicate"]), sep=":")))
  node_vector <- append(node_vector, as.character(unique(exp_design[,"Original File Name"])))
  node_vector <- append(node_vector, as.character(paste("Crux",unique(exp_design[,"Crux File Integer"]), sep=":")))
  node_shape <- rep("rectangular", length(node_vector))
  
  from_vector <- c()
  to_vector <- c()
  i <- 1
  for (i in 1:nrow(exp_design)){
    from_vector <- append(from_vector, as.character(exp_design[i,]$`Test or Control`))
    to_vector <- append(to_vector, as.character(exp_design[i,]$`Biological Condition`))
    
    from_vector <- append(from_vector, as.character(exp_design[i,]$`Biological Condition`))
    to_vector <- append(to_vector, paste("BioRep",as.character(exp_design[i,]$BioReplicate), sep=":"))
    
    from_vector <- append(from_vector, paste("BioRep",as.character(exp_design[i,]$BioReplicate), sep=":"))
    to_vector <- append(to_vector, as.character(exp_design[i,]$`Original File Name`))
    
    from_vector <- append(from_vector, as.character(exp_design[i,]$`Original File Name`))
    to_vector <- append(to_vector, paste("Crux", as.character(exp_design[i,]$`Crux File Integer`), sep=":"))
  }
  edge_df <- melt(data.frame(from_vector,to_vector))
  edge_df <- edge_df[!duplicated(edge_df[1:2]),]
  edge_df <- edge_df[order(edge_df["from_vector"]),]
  from_vector <- as.vector(edge_df[['from_vector']])
  to_vector <- as.vector(edge_df[['to_vector']])
  # Create a node data frame
  nodes <-
    create_nodes(nodes = node_vector,
                 shape = node_shape,
                 label = FALSE,
                 type = "lower",
                 color = "grey"
    )
  
  edges <-
    create_edges(from = from_vector,
                 to = to_vector,
                 relationship = "leading_to")
  
  
  graph <-
    create_graph(nodes_df = nodes,
                 edges_df = edges,
                 graph_attrs = "rankdir = LR",
                 node_attrs = "fontname = Helvetica",
                 edge_attrs = c("color = black",
                                "arrowsize = 2"))
  return(graph)
}