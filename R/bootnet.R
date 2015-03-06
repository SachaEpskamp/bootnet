### function to bootstrap networks ###
# Input: 
## data that can be bootstrapped
## a function to preprocess data into input for the network function
## A function that computes the network
## A function that extracts network, means and intercepts
## 
## The function will automatically detect if your data is continuous or binary and run 
## qgraph's EBICglasso or IsingFit appropriatly.

bootnet <- function(
  data, # Dataset
  nBoots = 1000, # Number of bootstrap samples.
  default, # Default method to use. EBICglasso, IsingFit, concentration, some more....
  prepFun, # Fun to produce the correlation or covariance matrix
  prepArgs = list(), # list with arguments for the correlation function
  estFun, # function that results in a network
  estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
  graphFun, # set to identity if missing
  graphArgs, # Set to null if missing
  intFun, # Set to null if missing
  intArgs, # Set to null if missing
  verbose = TRUE, # messages on what is being done?
  labels, # if missing taken from colnames
  alpha = 1# centrality alpha
){
  # First test if data is a data frame:
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  if (missing(labels)){
    labels <- colnames(data)
  }
  
  # Check and remove any variable that is not ordered, integer or numeric:
  goodColumns <- sapply(data, function(x) is.numeric(x) | is.ordered(x) | is.integer(x))
  
  if (!all(goodColumns)){
    if (verbose){
      warning(paste0("Removing non-numeric columns: ",paste(which(!goodColumns),collapse="; ")))
    }
    data <- data[,goodColumns,drop=FALSE]
  }
  
  
  ### DEFAULT OPTIONS ###
  if (missing(default)){
    if (missing(prepFun) | missing(prepArgs) | missing(estFun) | missing(estArgs)){
      stop("If 'default' is not set, 'prepFun', 'prepArgs', 'estFun' and 'estArgs' may not be missing.")
    }
  }

  if (!missing(default)){
    # prepFun:
    if (missing(prepFun)){
      prepFun <- switch(default,
                        EBICglasso = qgraph::cor_auto,
                        IsingFit = binarize,
                        pcor = qgraph::cor_auto
      )      
    }
    
    # prepArgs:
    qgraphVersion <- packageDescription("qgraph")$Version
    qgraphVersion <- as.numeric(strsplit(qgraphVersion,split="\\.|\\-")[[1]])
    if (length(qgraphVersion)==1) qgraphVersion <- c(qgraphVersion,0)
    if (length(qgraphVersion)==2) qgraphVersion <- c(qgraphVersion,0)
    goodVersion <- 
      (qgraphVersion[[1]] >= 1 & qgraphVersion[[2]] >= 3 & qgraphVersion[[3]] >= 1) | 
      (qgraphVersion[[1]] >= 1 & qgraphVersion[[2]] > 3) | 
      qgraphVersion[[1]] > 1
    
    newqgraph <- 
    if (missing(prepArgs)){
      prepArgs <- switch(default,
                         EBICglasso = ifelse(goodVersion,list(verbose=FALSE),list()),
                         IsingFit = list(),
                         pcor = list(verbose=FALSE)
      )
    }
    
    # estFun:
    if (missing(estFun)){
      estFun <- switch(default,
                       EBICglasso = qgraph::EBICglasso,
                       IsingFit = IsingFit::IsingFit,
                       pcor = corpcor::cor2pcor
      )
    }
    
    # estArgs:
    if (missing(estArgs)){
      estArgs <- switch(default,
                        EBICglasso = list(n = nrow(data), returnAllResults = TRUE),
                        IsingFit = list(plot = FALSE, progress = FALSE),
                        pcor = list()
      )
    }
    
    # graphFun:
    if (missing(graphFun)){
      graphFun <- switch(default,
                         EBICglasso = function(x)x[['optnet']],
                         IsingFit = function(x)x[['weiadj']],
                         pcor = identity
      )
    }
    
    # graphArgs:
    if (missing(graphArgs)){
      graphArgs <- switch(default,
                          EBICglasso = list(),
                          IsingFit = list(),
                          pcor = list()
      )
    }
    
    # intFun:
    if (missing(intFun)){
      intFun <- switch(default,
                       EBICglasso = null,
                       IsingFit = function(x)x[['thresholds']],
                       pcor = null
      )
    }
    
    
  }
  
  if (missing(prepFun)){
    prepFun <- identity
  }
  
  if (missing(prepArgs)){
    prepArgs <- list()
  }
  
  if (missing(graphFun)){
    graphFun <- identity
  }
  
  if (missing(graphArgs)){
    graphArgs <- list()
  }
  
  if (missing(intFun)){
    intFun <- null
  }
  
  if (missing(intArgs)){
    intArgs <- list()
  }
  
  # Estimate sample result:
  if (verbose){
    message("Estimating sample network...")
  }
  res <- estimateNetwork(data, prepFun, prepArgs, estFun, estArgs)
  sampleResult <- list(
    graph = do.call(graphFun,c(list(res), graphArgs)),
    intercepts = do.call(intFun,c(list(res), intArgs)),
    results = res,
    labels = labels
  )
  class(sampleResult) <- c("bootnetResult", "list")
  
  if (!isSymmetric(sampleResult[['graph']])){
    stop("bootnet does not support directed graphs")
  }
  
  # Bootstrap results:
  bootResults <- vector("list", nBoots)
  
  if (verbose){
    message("Bootstrapping...")
    pb <- txtProgressBar(0,nBoots,style = 3)
  }
  
  for (b in seq_len(nBoots)){
    bootData <- data[sample(seq_len(nrow(data)), nrow(data), replace=TRUE), ]
    res <- estimateNetwork(bootData, prepFun, prepArgs, estFun, estArgs)
    bootResults[[b]] <- list(
      graph = do.call(graphFun,c(list(res), graphArgs)),
      intercepts = do.call(intFun,c(list(res), intArgs)),
      results = res,
      labels = labels
    )
    class(bootResults[[b]]) <- c("bootnetResult", "list")
    
    if (verbose){
      setTxtProgressBar(pb, b)
    }
  }
  if (verbose){
    close(pb)
  }
  

  ### Compute the full parameter table!!
  if (verbose){
    message("Computing statistics...")
    pb <- txtProgressBar(0,nBoots+1,style = 3)
  }
  statTableOrig <- statTable(sampleResult,  name = "sample", alpha = alpha)
  if (verbose){
    setTxtProgressBar(pb, 1)
  }
  statTableBoots <- vector("list", nBoots)
  for (b in seq_len(nBoots)){
    statTableBoots[[b]] <- statTable(bootResults[[b]], name = paste("boot",b), alpha = alpha)    
    if (verbose){
      setTxtProgressBar(pb, b+1)
    }
  }
  if (verbose){
    close(pb)
  }
  
  # Ordereing by node name to make nice paths:

  Result <- list(
    sampleTable = statTableOrig,
    bootTable =  dplyr::rbind_all(statTableBoots),
    sample = sampleResult,
    boots = bootResults)
  
  class(Result) <- "bootnet"
  
  return(Result)
}