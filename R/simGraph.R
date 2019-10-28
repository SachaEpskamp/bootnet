genGGM <- function(
  Nvar,
  p = 0, # Rewiring probability if graph = "smallworld" or "cluster", or connection probability if graph = "random". If cluster, can add multiple p's for each cluster, e.g., "c(.1, .5)" 
  nei = 1,
  parRange = c(0.5,1),
  constant = 1.5,
  propPositive = 0.5,
  clusters = NULL, #number of clusters if graph = "cluster"
  graph = c("smallworld","random", "scalefree", "hub", "cluster")
){
  graph <- match.arg(graph)
  
  
  ## Approach from 
  # Yin, J., & Li, H. (2011). A sparse conditional gaussian graphical model for analysis of genetical genomics data. The annals of applied statistics, 5(4), 2630.
  
  # Simulate graph structure:
  if (graph == "smallworld"){
    # Watts Strogatz small-world
    trueKappa <- as.matrix(igraph::get.adjacency(igraph::watts.strogatz.game(1,Nvar,nei,p)))    
  } else if (graph == "random"){
    # Ranodm network:
    trueKappa <- as.matrix(igraph::get.adjacency(igraph::erdos.renyi.game(Nvar, p)))
  } else if (graph == "scalefree") {
    trueKappa <- BDgraph::bdgraph.sim(p = Nvar, graph = "scale-free")$G
  } else if (graph == "hub") {
    trueKappa <- BDgraph::bdgraph.sim(p = Nvar, graph = "hub")$G
    class(trueKappa) <- "matrix"
  } else if (graph == "cluster") {
    trueKappa <-  BDgraph::bdgraph.sim(p = Nvar, graph = "cluster", prob = p, class = clusters)$G #can be 
    class(trueKappa) <- "matrix"
  }
  
  # Make edges negative and add weights:
  trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] * sample(c(-1,1),sum(upper.tri(trueKappa)),TRUE,prob=c(propPositive,1-propPositive)) * 
    runif(sum(upper.tri(trueKappa)), min(parRange ),max(parRange ))
  
  # Symmetrize:
  trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]  
  
  # Make pos def:
  diag(trueKappa) <- constant * rowSums(abs(trueKappa))
  diag(trueKappa) <- ifelse(diag(trueKappa)==0,1,diag(trueKappa))
  trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
  trueKappa <- (trueKappa + t(trueKappa)) / 2
  
  return(as.matrix(qgraph::wi2net(trueKappa)))
}

# genGGM <- function(
#   Nvar,
#   p = 0, # Rewiring probability if graph = "smallworld", or connection probability if graph = "random"
#   nei = 1,
#   parRange = c(0.5,1),
#   constant = 1.5,
#   propPositive = 0.5,
#   graph = c("smallworld","random")
# ){
#   graph <- match.arg(graph)
# 
#   
#   ## Approach from 
#   # Yin, J., & Li, H. (2011). A sparse conditional gaussian graphical model for analysis of genetical genomics data. The annals of applied statistics, 5(4), 2630.
#   
#   # Simulate graph structure:
#   if (graph == "smallworld"){
#     # Watts Strogatz small-world
#     trueKappa <- as.matrix(igraph::get.adjacency(igraph::watts.strogatz.game(1,Nvar,nei,p)))    
#   } else if (graph == "random"){
#     # Ranodm network:
#     trueKappa <- as.matrix(igraph::get.adjacency(igraph::erdos.renyi.game(Nvar, p)))
#   }
# 
#   
#   # Make edges negative and add weights:
#   trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] * sample(c(-1,1),sum(upper.tri(trueKappa)),TRUE,prob=c(propPositive,1-propPositive)) * 
#     runif(sum(upper.tri(trueKappa)), min(parRange ),max(parRange ))
#   
#   # Symmetrize:
#   trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]  
#   
#   # Make pos def:
#   diag(trueKappa) <- constant * rowSums(abs(trueKappa))
#   diag(trueKappa) <- ifelse(diag(trueKappa)==0,1,diag(trueKappa))
#   trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
#   trueKappa <- (trueKappa + t(trueKappa)) / 2
#   
#   return(as.matrix(qgraph::wi2net(trueKappa)))
# }