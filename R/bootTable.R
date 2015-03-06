# Compute tidy table from bootnetResult object:
# Result in data frame with entries:
# original (logical)
# name
# type
# node1
# node2
# value

statTable <- function(x, name, alpha = 1){
  stopifnot(is(x, "bootnetResult"))
  tables <- list()
  
  # edges:
  ind <- which(upper.tri(x[['graph']], diag=FALSE), arr.ind=TRUE)
  tables$edges <- dplyr::tbl_df(data.frame(
    name = name,
    type = "edge",
    node1 = x[['labels']][ind[,1]],
    node2 = x[['labels']][ind[,2]],
    value = x[['graph']][upper.tri(x[['graph']], diag=FALSE)]
    ))
  
  # Intercepts:
  if (!is.null(x[['intercepts']])){
    tables$intercepts <- dplyr::tbl_df(data.frame(
      name = name,
      type = "intercept",
      node1 = x[['labels']],
      node2 = '',
      value = x[['intercepts']]
    ))
  } 
  # Centrality analysis:
  cent <- qgraph::centrality(x[['graph']], alpha = alpha)

  # strength:
  tables$strength <- dplyr::tbl_df(data.frame(
    name = name,
    type = "strength",
    node1 = x[['labels']],
    node2 = '',
    value = cent[['OutDegree']]
  ))
  
  # closeness:
  tables$closeness <- dplyr::tbl_df(data.frame(
    name = name,
    type = "closeness",
    node1 = x[['labels']],
    node2 = '',
    value = cent[['Closeness']]
  ))
  
  
  # betweenness:
  tables$betweenness <- dplyr::tbl_df(data.frame(
    name = name,
    type = "betweenness",
    node1 = x[['labels']],
    node2 = '',
    value = cent[['Betweenness']]
  ))
  
  tables$sp <- dplyr::tbl_df(data.frame(
    name = name,
    type = "distance",
    node1 = x[['labels']][ind[,1]],
    node2 = x[['labels']][ind[,2]],
    value = cent[['ShortestPathLengths']][upper.tri(cent[['ShortestPathLengths']], diag=FALSE)]
  ))
  
  for (i in seq_along(tables)){
    tables[[i]]$id <- ifelse(tables[[i]]$node2=='',paste0("N: ",tables[[i]]$node1),paste0("E: ",tables[[i]]$node1, "--", tables[[i]]$node2))
  }
  
  return(dplyr::rbind_all(tables))
}
