genGGM <- function(
  Nvar,
  p = 0,
  nei = 1,
  parRange = c(0.5,1),
  constant = 1.5,
  propPositive = 0.5
){
  # Watts Strogatz small-world
  
  ## Approach from 
  # Yin, J., & Li, H. (2011). A sparse conditional gaussian graphical model for analysis of genetical genomics data. The annals of applied statistics, 5(4), 2630.
  
  # Empty matrix:
  #   trueKappa <- matrix(0,Nvar,Nvar)
  #   
  #   # Total edges:
  #   totEdges <- sum(upper.tri(trueKappa))
  #   
  #   # Included edges:
  #   nEdges <- round((1-sparsity)*totEdges)
  #   
  #   # Sample the edges:
  #   inclEdges <- sample(seq_len(totEdges),nEdges)
  #   
  #   # Make edges:
  #   trueKappa[upper.tri(trueKappa)][inclEdges] <- 1
  trueKappa <- as.matrix(igraph::get.adjacency(igraph::watts.strogatz.game(1,Nvar,nei,p)))
  
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