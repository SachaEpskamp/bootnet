# Function to run a standard network simulation study:
# Two methods: dataGenerator function, or using graph, intercepts and model for GGM/Ising.

# Generator functions:
ggmGenerator <- function(
  ordinal = FALSE,
  nLevels = 4
){
  # Generate a function with nCase as first and input as second argument
  
  # Wrapper:
  Estimator <- function(n, input){
    if (is.list(input)){
      if ("graph" %in% names(input)){
        graph <- input$graph
      } else stop("'graph' not in input list.")
      
      if ("intercepts" %in% names(input)){
        intercepts <- input$intercepts
      } else {
        intercepts <- rep(0,ncol(graph))
      }
    } else {
      if (!is.matrix(input)){
        stop("'input' is not a matrix or list.")
      }
      
      graph <- input
      intercepts <- rep(0,ncol(graph))
    }
    # standardize:
    if (!all(diag(graph) == 0 | diag(graph) == 1)){
      graph <- cov2cor(graph)
    }
    
    # Remove diag:
    diag(graph) <- 0
    
    # Generate data:
    # True sigma:
    Sigma <- cov2cor(solve(diag(ncol(graph)) - graph))
    
    # Generate data:
    Data <- mvtnorm::rmvnorm(n, sigma = Sigma)
    
    ORDINALDUMMY
    
    return(Data)
  }
  
  # Deparse:
  deparsedEstimator <- deparse(Estimator)
  l <- grep("ORDINALDUMMY",deparsedEstimator)
  
  if (!ordinal){
    deparsedEstimator <- deparsedEstimator[-l]
  } else {
    deparsedEstimator[l] <- sprintf('
    for (i in 1:ncol(Data)){
      Data[,i] <- as.numeric(cut(Data[,i],sort(c(-Inf,rnorm(%f-1),Inf))))
    }',nLevels)
  }

  # Parse again:
  Estimator <- eval(parse(text=deparsedEstimator))
  return(Estimator)
}


netSimulator <- function(
  input = simGraph(Nvar = 10),  # A matrix, or a list with graph and intercepts elements.
  dataGenerator = ggmGenerator(ordinal = FALSE), # Data generating function taking N as first argument and input as second.
  nCases = c(50,100,250,500,1000,2500), # Number of cases
  nReps = 100, # Number of repititions per condition
  nCores = 1, # Number of computer cores used
  ... # estimateNetwork arguments (if none specified, will default to default = "EBICglasso)
){
  # Dots list:
  .dots <- list(...)
  if (length(.dots) == 0){
    message("Setting default = 'EBICglasso'")
    .dots$default <- "EBICglasso"
  }

  Results <- parSim(
    # Conditions:
    nCases = nCases,
    
    # Setup:,
    write=FALSE,
    nCores = nCores,
    reps = nReps,
    debug=FALSE,
    export = c("input","dataGenerator",".dots"),
    
    expression = {
      cor0 <- function(x,y){
        if (all(is.na(x)) || all(is.na(y))){
          return(NA)
        }
        
        if (sd(x,na.rm=TRUE)==0 | sd(y,na.rm=TRUE) == 0){
          return(0)
        }
        
        return(cor(x,y,use="pairwise.complete.obs"))
      }
      
      # Generate the input:
      if (is.function(input)){
        inputResults <- input()
      } else {
        inputResults <- input
      }
      
      # True network:
      if (is.list(inputResults)){
        trueNet <- inputResults$graph
      } else {
        trueNet <- inputResults
      }
      
      # Generate the data:
      Data <- dataGenerator(nCases, inputResults)

      # Compute the network:
      .dots$data <- Data
      netResults <- do.call(bootnet::estimateNetwork,.dots)
      estNet <- qgraph::getWmat(netResults)
      
      # Compute measures:
      ### STORE RESULTS ###
      SimulationResults <- list()
      
      # Estimated edges:
      est <- estNet[upper.tri(estNet)]
      # Real edges:
      real <- trueNet[upper.tri(trueNet)]
      
      # Equal?
      SimulationResults$correctModel <- all((est == 0) == (real == 0))
      
      # True positives:
      TruePos <- sum(est != 0 &  real != 0)
      
      # False pos:
      FalsePos <- sum(est != 0 & real == 0)
      
      # True Neg:
      TrueNeg <- sum(est == 0 & real == 0)
      
      # False Neg:
      FalseNeg <- sum(est == 0 & real != 0)
      
      ### Sensitivity:
      SimulationResults$sensitivity <- TruePos / (TruePos + FalseNeg)
      
      # Specificity:
      SimulationResults$specificity <- TrueNeg / (TrueNeg + FalsePos)
      
      # Correlation:
      SimulationResults$correlation <- cor0(est,real)
      
      SimulationResults
    })
  
  class(Results) <- c("netSimulator","data.frame")
  return(Results)
}

