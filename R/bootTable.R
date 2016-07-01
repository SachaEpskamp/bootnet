# Compute tidy table from bootnetResult object:
# Result in data frame with entries:
# original (logical)
# name
# type
# node1
# node2
# value

statTable <- function(x, name, alpha = 1, computeCentrality = TRUE){
  type <- NULL
  value <- NULL
  
  stopifnot(is(x, "bootnetResult"))
  tables <- list()
  
  # edges:
  ind <- which(upper.tri(x[['graph']], diag=FALSE), arr.ind=TRUE)
  tables$edges <- dplyr::tbl_df(data.frame(
    name = name,
    type = "edge",
    node1 = x[['labels']][ind[,1]],
    node2 = x[['labels']][ind[,2]],
    value = x[['graph']][upper.tri(x[['graph']], diag=FALSE)],
    stringsAsFactors = FALSE
    ))
  
  
  tables$length <- dplyr::tbl_df(data.frame(
    name = name,
    type = "length",
    node1 = x[['labels']][ind[,1]],
    node2 = x[['labels']][ind[,2]],
    value = abs(1/abs(x[['graph']][upper.tri(x[['graph']], diag=FALSE)])),
    stringsAsFactors = FALSE
  ))
  
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
    cent <- qgraph::centrality(x[['graph']], alpha = alpha)
    
    # strength:
    tables$strength <- dplyr::tbl_df(data.frame(
      name = name,
      type = "strength",
      node1 = x[['labels']],
      node2 = '',
      value = cent[['OutDegree']],
      stringsAsFactors = FALSE
    ))
    
    # closeness:
    tables$closeness <- dplyr::tbl_df(data.frame(
      name = name,
      type = "closeness",
      node1 = x[['labels']],
      node2 = '',
      value = cent[['Closeness']],
      stringsAsFactors = FALSE
    ))
    
    
    # betweenness:
    tables$betweenness <- dplyr::tbl_df(data.frame(
      name = name,
      type = "betweenness",
      node1 = x[['labels']],
      node2 = '',
      value = cent[['Betweenness']],
      stringsAsFactors = FALSE
    ))
    
    tables$sp <- dplyr::tbl_df(data.frame(
      name = name,
      type = "distance",
      node1 = x[['labels']][ind[,1]],
      node2 = x[['labels']][ind[,2]],
      value = cent[['ShortestPathLengths']][upper.tri(cent[['ShortestPathLengths']], diag=FALSE)],
      stringsAsFactors = FALSE
    ))
    
  
  }
  for (i in seq_along(tables)){
    tables[[i]]$id <- ifelse(tables[[i]]$node2=='',paste0("N: ",tables[[i]]$node1),paste0("E: ",tables[[i]]$node1, "--", tables[[i]]$node2))
  }  
  
  tab <- dplyr::rbind_all(tables)
  tab$nNode <- x$nNodes
  tab$nPerson <- x$nPerson
  
  # Compute rank:
  tab <- tab %>% group_by(type) %>%
    mutate(rank_avg = rank(value,ties.method = "average"),
           rank_min = rank(value,ties.method = "min"),
           rank_max = rank(value,ties.method = "max"))
  
  return(tab)
}
