# Compute tidy table from bootnetResult object:
# Result in data frame with entries:
# original (logical)
# name
# type
# node1
# node2
# value

statTable <- function(x, name, alpha = 1, computeCentrality = TRUE,statistics = c("edge","strength","closeness","betweenness")){
  # Statistics can be:
  if (!all(statistics %in% c("intercept","edge","length","distance","closeness","betweenness","strength"))){
    stop("'statistics' must be 'edge', 'intercept', 'length', 'distance', 'closeness', 'betweenness' or 'strength'")
  }
  
  
  type <- NULL
  value <- NULL
  
  stopifnot(is(x, "bootnetResult"))
  tables <- list()
  if (is.null(x[['labels']])){
    x[['labels']] <- seq_len(ncol(x[['graph']]))
  }
  
  # edges:
  ind <- which(upper.tri(x[['graph']], diag=FALSE), arr.ind=TRUE)
  
  # Weights matrix:
  Wmat <- qgraph::getWmat(x)
  
  if ("edge" %in% statistics){
    tables$edges <- dplyr::tbl_df(data.frame(
      name = name,
      type = "edge",
      node1 = x[['labels']][ind[,1]],
      node2 = x[['labels']][ind[,2]],
      value = Wmat[upper.tri(Wmat, diag=FALSE)],
      stringsAsFactors = FALSE
    ))
  }
  
  
  if ("length" %in% statistics){ 
    tables$length <- dplyr::tbl_df(data.frame(
      name = name,
      type = "length",
      node1 = x[['labels']][ind[,1]],
      node2 = x[['labels']][ind[,2]],
      value = abs(1/abs(Wmat[upper.tri(Wmat, diag=FALSE)])),
      stringsAsFactors = FALSE
    ))
  }
  
  # Intercepts:
  if (!is.null(x[['intercepts']])){
    tables$intercepts <- dplyr::tbl_df(data.frame(
      name = name,
      type = "intercept",
      node1 = x[['labels']],
      node2 = '',
      value = x[['intercepts']],
      stringsAsFactors = FALSE
    ))
  } 
  
  if (computeCentrality){
    # Centrality analysis:
    if (all(x[['graph']]==0)){
      cent <- list(
        OutDegree = rep(0,ncol(x[['graph']])),
        InDegree = rep(0,ncol(x[['graph']])),
        Closeness = rep(0,ncol(x[['graph']])),
        Betweenness = rep(0,ncol(x[['graph']])),
        ShortestPathLengths = matrix(Inf,ncol(x[['graph']]),ncol(x[['graph']]))
      )
    } else {
      cent <- qgraph::centrality(Wmat, alpha = alpha, all.shortest.paths = FALSE)
      
    }
    
    # strength:
    if ("strength" %in% statistics){
      
    tables$strength <- dplyr::tbl_df(data.frame(
      name = name,
      type = "strength",
      node1 = x[['labels']],
      node2 = '',
      value = cent[['OutDegree']],
      stringsAsFactors = FALSE
    ))
    }
    
    # closeness:
      if ("closeness" %in% statistics){
    tables$closeness <- dplyr::tbl_df(data.frame(
      name = name,
      type = "closeness",
      node1 = x[['labels']],
      node2 = '',
      value = cent[['Closeness']],
      stringsAsFactors = FALSE
    ))
      }
    
    
    # betweenness:
    if ("betweenness" %in% statistics){
    tables$betweenness <- dplyr::tbl_df(data.frame(
      name = name,
      type = "betweenness",
      node1 = x[['labels']],
      node2 = '',
      value = cent[['Betweenness']],
      stringsAsFactors = FALSE
    ))
    }
    
    if ("distance" %in% statistics){
    tables$sp <- dplyr::tbl_df(data.frame(
      name = name,
      type = "distance",
      node1 = x[['labels']][ind[,1]],
      node2 = x[['labels']][ind[,2]],
      value = cent[['ShortestPathLengths']][upper.tri(cent[['ShortestPathLengths']], diag=FALSE)],
      stringsAsFactors = FALSE
    ))
    }
    
    
  }
  #   for (i in seq_along(tables)){
  #     tables[[i]]$id <- ifelse(tables[[i]]$node2=='',paste0("N: ",tables[[i]]$node1),paste0("E: ",tables[[i]]$node1, "--", tables[[i]]$node2))
  #   }  
  
  for (i in seq_along(tables)){
    tables[[i]]$id <- ifelse(tables[[i]]$node2=='',tables[[i]]$node1,paste0(tables[[i]]$node1, "--", tables[[i]]$node2))
  }  
  
  tab <- dplyr::bind_rows(tables)
  tab$nNode <- x$nNodes
  tab$nPerson <- x$nPerson
  
  # Compute rank:
  tab <- tab %>% group_by(type) %>%
    mutate(rank_avg = rank(value,ties.method = "average"),
           rank_min = rank(value,ties.method = "min"),
           rank_max = rank(value,ties.method = "max"))
  
  return(tab)
}
