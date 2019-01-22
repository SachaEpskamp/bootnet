# Function to create include probability network
bootInclude <- function(bootobject,verbose=TRUE){
  # Check if object is bootnet object:
  if (class(bootobject) != "bootnet"){
    stop("'bootobject' must be an object of class 'bootnet'")
  }
  
  # Check type:
  if (bootobject$type != "nonparametric" & bootobject$type != "parametric"){
    stop("Bootstrap type must be 'nonparametric' or 'parametric'")
  }
  
  # Extract the network object:
  Network <- bootobject$sample
  # Dummy for multiple graphs:
  if (!is.list(Network$graph)){
    Graphs <- list(Network$graph)
    Directed <- list(Network$directed)
    Intercepts <- list(Network$intercepts)
    names(Graphs) <- names(Directed) <- names(Intercepts) <-
      unique(bootobject$bootTable$graph)
  } else {
    Graphs <- Network$graph
    Directed <- Network$directed
    Intercepts <-  Network$intercepts
  }
  
  # For every graph:
  for (g in seq_along(Graphs)){
    graphName <- names(Graphs)[g]

    # Summary table of edge weights:
    bootSummary <- bootobject$bootTable %>% 
      dplyr::filter_(~type == "edge", ~graph == graphName) %>%
      dplyr::group_by_(~node1,~node2) %>%
      dplyr::summarize(
        propNonZero=mean(value != 0)
      )
    
    # Reweight network:
    # if (nrow(bootSummary) > 0){
    Graphs[[graphName]][] <- 0

    
    for (i in 1:nrow(bootSummary)){
      Graphs[[graphName]][Network$labels == bootSummary$node1[i],Network$labels == bootSummary$node2[i]] <- bootSummary$propNonZero[i]
      if (!Directed[[graphName]]){
        Graphs[[graphName]][Network$labels == bootSummary$node2[i],Network$labels == bootSummary$node1[i]] <- bootSummary$propNonZero[i]
      }
    }
  }
  
  # Return to network object:
  if (length(Graphs) == 1){
    Network$graph <- Graphs[[1]]
    Network$intercepts <- NULL
  } else {
    Network$graph <- Graphs
    Network$intercepts <- NULL
  }
  
  # Add indicator network is about include proportions:
  Network$bootInclude <- TRUE
  
  # Return network:
  return(Network)
}