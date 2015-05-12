# This function takes data as input and produced a network. It is used inside bootnet:
estimateNetwork <- function(
  data,
  prepFun, # Fun to produce the correlation or covariance matrix
  prepArgs = list(), # list with arguments for the correlation function
  estFun, # function that results in a network
  estArgs # arguments sent to the graph estimation function (if missing automatically sample size is included)
){
  # Compute input:
  input <- do.call(prepFun, c(list(data), prepArgs))
  
  # Compute network:
  Res <- do.call(estFun, c(list(input),estArgs))
  
  # Return network:
  return(Res)
}