# Function to run a standard network simulation study:
# Two methods: dataGenerator function, or using graph, intercepts and model for GGM/Ising.

# Generator functions:
ggmGenerator <- function(
  ordinal = FALSE,
  nLevels = 4
){
  ORDINALDUMMY <- NULL
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
    if (any(eigen(diag(ncol(graph)) - graph)$values < 0)){
      stop("Precision matrix is not positive semi-definite")
    }
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
         if (is.list(input) && !is.null(input$thresholds)) {
        Data[,i] <- as.numeric(cut(Data[,i],sort(c(-Inf,input$thresholds[[i]],Inf))))
    } else {
        Data[,i] <- as.numeric(cut(Data[,i],sort(c(-Inf,rnorm(%f-1),Inf))))
    }  
  
    }',nLevels)
  }
  
  # Parse again:
  Estimator <- eval(parse(text=deparsedEstimator))
  return(Estimator)
}

### Ising:
IsingGenerator <- function(
  ... # Arguments used in IsingSampler
){
  ARGDUMMY <- NULL
  # Generate a function with nCase as first and input as second argument
  
  # Wrapper:
  Gen <- function(n, input){
    if (is.list(input)){
      if ("graph" %in% names(input)){
        graph <- input$graph
      } else stop("Input must be a list containing elements 'graph' and 'intercepts'")
      
      if ("intercepts" %in% names(input)){
        intercepts <- input$intercepts
      } else {
        stop("Input must be a list containing elements 'graph' and 'intercepts'")
      }
    } else {
      stop("Input must be a list containing elements 'graph' and 'intercepts'")
    }
    
    Data <- IsingSampler::IsingSampler(n, 
                                       input$graph, 
                                       ARGDUMMY,
                                       input$intercepts
                                       
    )
    
    return(Data)
  }
  
  # Deparse:
  deparsedEstimator <- deparse(Gen)
  l <- grep("ARGDUMMY,",deparsedEstimator)
  
  dots <-list(...)
  
  
  
  if (length(dots) == 0){
    deparsedEstimator[l] <- gsub("ARGDUMMY,","",deparsedEstimator[l])
  } else {
    # Construct arguments:
    txt <- sapply(dots, function(x)paste(deparse(dput(x)),collapse="\n"))
    deparsedEstimator[l] <- gsub("ARGDUMMY",paste(names(dots), "=", txt, collapse = ", "),deparsedEstimator[l])
  }
  
  # Parse again:
  Estimator <- eval(parse(text=deparsedEstimator))
  return(Estimator)
}


netSimulator <- function(
  input = genGGM(Nvar = 10),  # A matrix, or a list with graph and intercepts elements. Or a generating function
  nCases = c(50,100,250,500,1000,2500), # Number of cases
  nReps = 100, # Number of repititions per condition
  nCores = 1, # Number of computer cores used
  default,
  dataGenerator, 
  ..., # estimateNetwork arguments (if none specified, will default to default = "EBICglasso)
  moreArgs = list(), # List of extra args not intended to be varied as conditions
  moreOutput = list() # List with functions that take two weights matrices and produce some value
){
  # Dots list:
  .dots <- list(...)
  # default <- match.arg(default)
  
  # Check default and dataGenerator:
  if (missing(default) & missing(dataGenerator)){
    message("'default' and 'dataGenerator' are missing. Setting default = 'EBICglasso'")
    default <- "EBICglasso"
  }
  
  if (missing(default) & length(.dots) == 0){
    message("No estimator specified. Setting default = 'EBICglasso'")
    default <- "EBICglasso"
  }
  
  # Data generator:
  if (!missing(default) & missing(dataGenerator)){
    
    if (default == "EBICglasso" || default == "glasso" || default == "pcor" || default == "adalasso" || default == "huge"|| default == "ggmModSelect" || default == "LoGo"){
      message("Setting 'dataGenerator = ggmGenerator(ordinal = FALSE)'")
      dataGenerator <- ggmGenerator(ordinal = FALSE)
    } else if (default == "IsingFit" || default == "IsingSampler"){
      message("Setting 'dataGenerator = IsingGenerator()'")
      dataGenerator <- IsingGenerator()
    } else {
      stop(paste0("Default set '",default, "' not yet supported. Please manually specify 'dataGenerator'"))
    }
    
    
  }
  
  
  # Else none:
  if (missing(default)) default <- "none"
  
  # parSim arguments:
  Args <- c(
    list(
      # Conditions:
      nCases = nCases,
      default = default,
      
      # Setup:,
      write=FALSE,
      nCores = nCores,
      reps = nReps,
      debug=FALSE,
      export = c("input","dataGenerator",".dots","moreArgs","moreOutput"),
      
      expression = expression({
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
        args <- list()
        args$data <- Data
        args$verbose <- FALSE
        args$default <- default
        args <- c(args,moreArgs)
        
        for (i in seq_along(.dots)){
          args[[names(.dots)[i]]] <- get(names(.dots)[i])
        }
        netResults <- do.call(bootnet::estimateNetwork,args)
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
        
        # Centrality:
        centTrue <- qgraph::centrality(trueNet)
        centEst <- qgraph::centrality(estNet)
        
        SimulationResults$strength <- cor0(centTrue$OutDegree,centEst$OutDegree)
        SimulationResults$closeness <- cor0(centTrue$Closeness,centEst$Closeness)
        SimulationResults$betweenness <- cor0(centTrue$Betweenness,centEst$Betweenness)
        SimulationResults$ExpectedInfluence <- cor0(centTrue$OutExpectedInfluence,centEst$OutExpectedInfluence)

        # 
        # ### TEMP: REMOVE:
        # SimulationResults$MeanBiasFalsePositives <- mean(abs(est[real==0 & est!=0]))
        # # SimulationResults$Q75BiasFalsePositives <- quantile(abs(est[real==0 & est!=0]), 0.75)
        # SimulationResults$MaxBiasFalsePositives <- max(abs(est[real==0 & est!=0]))
        # SimulationResults$MaxWeight <- max(abs(est))
        # SimulationResults$MeanWeight <- mean(abs(est[est!=0]))
        if (any(real==0 & est!=0)){
          SimulationResults$MaxFalseEdgeWidth <- max(abs(est[real==0 & est!=0])) / max(abs(est))
        } else {
          SimulationResults$MaxFalseEdgeWidth <- NA
        }

        SimulationResults$bias <- mean(abs(est - real))
        # ##
        if (length(moreOutput) > 1){
          if (is.null(names(moreOutput))){
            names(moreOutput) <- paste0("moreOutput",seq_along(moreOutput))
          }
          
          for (out in seq_along(moreOutput)){
            SimulationResults[[out]] <- moreOutput[[i]](estNet, trueNet)
          }
        }
        
        SimulationResults
      })),
      .dots)

    Results <- do.call(parSim, Args)
    
    class(Results) <- c("netSimulator","data.frame")
    return(Results)
}

