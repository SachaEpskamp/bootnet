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
  tables$edges <- dplyr::data_frame(
    name = name,
    type = "edge",
    node1 = x[['labels']][ind[,1]],
    node2 = x[['labels']][ind[,2]],
    value = x[['graph']][upper.tri(x[['graph']], diag=FALSE)]
    )
  
  # Intercepts:
  if (!is.null(x[['intercepts']])){
    tables$intercepts <- dplyr::data_frame(
      name = name,
      type = "intercept",
      node1 = x[['labels']],
      node2 = '',
      value = x[['intercepts']]
    )
  } 
  # Centrality analysis:
  cent <- qgraph::centrality(x[['graph']], alpha = alpha)

  # strength:
  tables$strength <- dplyr::data_frame(
    name = name,
    type = "strength",
    node1 = x[['labels']],
    node2 = '',
    value = cent[['OutDegree']]
  )
  
  # closeness:
  tables$closeness <- dplyr::data_frame(
    name = name,
    type = "closeness",
    node1 = x[['labels']],
    node2 = '',
    value = cent[['Closeness']]
  )
  
  
  # betweenness:
  tables$betweenness <- dplyr::data_frame(
    name = name,
    type = "betweenness",
    node1 = x[['labels']],
    node2 = '',
    value = cent[['Betweenness']]
  )
  
  tables$sp <- dplyr::data_frame(
    name = name,
    type = "distance",
    node1 = x[['labels']][ind[,1]],
    node2 = x[['labels']][ind[,2]],
    value = cent[['ShortestPathLengths']][upper.tri(cent[['ShortestPathLengths']], diag=FALSE)]
  )
  
  for (i in seq_along(tables)){
    tables[[i]]$id <- ifelse(tables[[i]]$node2=='',paste0("N: ",tables[[i]]$node1),paste0("E: ",tables[[i]]$node1, "--", tables[[i]]$node2))
  }
  
  return(dplyr::rbind_all(tables))
}



# Wrong code:

# # Inner function to create a table from a vector of bootstraps of some parameter:
# bootnetTableEntry <- function(type, node1, node2, original, bootstraps, percentiles = c(1,2.5,5,25,50,75,95,97.5,99)){
#   if (length(original)!=1) stop("Original parameter must have length=1")
#   
#   # Compute boot quantiles:
#   bootPercentiles <- quantile(bootstraps, probs = percentiles/100)
#   df <- cbind(data.frame(type=type,node1=node1,node2=node2,orig=original,var=var(bootstraps),prop0 = mean(bootstraps==0)),
#               as.data.frame(t(bootPercentiles)))
#   rownames(df) <- NULL
#   return(df)
# }
# 
# # bootnetTable:
# bootnetTable <- function(
#   x, #bootnet result
#   type = c("edges","intercepts","centrality"),
#   percentiles = c(1,2.5,5,25,50,75,95,97.5,99),
#   diag = FALSE){
#   
#   stopifnot(is(x, "bootnet"))
#   # Gather statistics in a list of bootstrap results and original:
#   statList <- list()
# 
#   # Edges:
#   if ("edges" %in% type){
#     ind <- which(upper.tri(x[['sample']][['graph']], diag=diag), arr.ind=TRUE)
#     c <- length(statList)
#     for (i in seq_len(nrow(ind))){
#       statList[[c+i]] <- list(
#         type = "edge",
#         node1 = x[['labels']][ind[i,1]],
#         node2 = x[['labels']][ind[i,2]],
#         original = x[['sample']][['graph']][ind[i,1],ind[i,2]],
#         bootstraps = sapply(x[['boots']], function(x)x[['graph']][ind[i,1],ind[i,2]])
#       )
#     }
#   }
#   
#   # Intercepts:
#   if ("intercept" %in% type){
#     if (!is.null(x[['sample']][['intercept']])){
#       nInt <- length(x[['sample']][['intercept']])
#       c <- length(statList)
#       for (i in nInt){
#         statList[[c+i]] <- list(
#           type = "intercept",
#           node1 = x[['labels']][i],
#           node2 = '',
#           original = x[['sample']][['intercepts']][i],
#           bootstraps = sapply(x[['boots']], function(x)x[['intercepts']][i])
#         )
#       }
#       
#     }
#     
#     
#   }
# 
#   # construct table:
#   Table <- do.call(rbind, lapply(statList,function(x)do.call(bootnetTableEntry, c(x, list(percentiles = percentiles)))))
#   browser()
#   
# }