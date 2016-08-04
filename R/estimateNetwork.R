# This function takes data as input and produced a network. It is used inside bootnet:
estimateNetwork <- function(
  data,
  default = c("none", "EBICglasso", "pcor","IsingFit","IsingLL", "huge","adalasso"),
  prepFun, # Fun to produce the correlation or covariance matrix
  prepArgs = list(), # list with arguments for the correlation function
  estFun, # function that results in a network
  estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
  graphFun, # set to identity if missing
  graphArgs, # Set to null if missing
  intFun, # Set to null if missing
  intArgs, # Set to null if missing
  labels # if missing taken from colnames
  # plot = TRUE, # Plot the network?
  # ... # qgraph arguments
){
#   if (default[[1]]=="glasso") default <- "EBICglasso"
#   default <- match.arg(default)
#   
  # If NAs and default can't handle, stop:
  if (any(is.na(data)) && default %in% c("huge","adalasso")){
    stop(paste0("Missing data not supported for default set '",default,"'. Try using na.omit(data)."))
  }
  

  N <- ncol(data)
  Np <- nrow(data)
  
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
  
#   ### DEFAULT OPTIONS ###
#   if ((default == "none")){
#     if (missing(prepFun) | missing(prepArgs) | missing(estFun) | missing(estArgs)){
#       stop("If 'default' is not set, 'prepFun', 'prepArgs', 'estFun' and 'estArgs' may not be missing.")
#     }
#   }
#   
#   if (!(default == "none")){
#     # prepFun:
#     if (missing(prepFun)){
#       prepFun <- switch(default,
#                         EBICglasso = qgraph::cor_auto,
#                         IsingFit = binarize,
#                         IsingLL = binarize,
#                         pcor = qgraph::cor_auto,
#                         huge = function(x)huge::huge.npn(na.omit(as.matrix(x)),verbose = FALSE),
#                         adalasso = identity
#       )
#       #       prepFun <- switch(default,
#       #                         EBICglasso = cor,
#       #                         IsingFit = binarize,
#       #                         pcor = cor
#       #       )      
#     }
#     
#     # prepArgs:
#     #     qgraphVersion <- packageDescription("qgraph")$Version
#     #     qgraphVersion <- as.numeric(strsplit(qgraphVersion,split="\\.|\\-")[[1]])
#     #     if (length(qgraphVersion)==1) qgraphVersion <- c(qgraphVersion,0)
#     #     if (length(qgraphVersion)==2) qgraphVersion <- c(qgraphVersion,0)
#     #     goodVersion <- 
#     #       (qgraphVersion[[1]] >= 1 & qgraphVersion[[2]] >= 3 & qgraphVersion[[3]] >= 1) | 
#     #       (qgraphVersion[[1]] >= 1 & qgraphVersion[[2]] > 3) | 
#     #       qgraphVersion[[1]] > 1
#     
#     
#     if (missing(prepArgs)){
#       prepArgs <- switch(default,
#                          EBICglasso = ifElse(identical(prepFun,qgraph::cor_auto),list(verbose=FALSE),
#                                              ifElse(identical(prepFun,cor),list(use = "pairwise.complete.obs"),list())),
#                          IsingFit = list(),
#                          pcor =  ifElse(identical(prepFun,qgraph::cor_auto),list(verbose=FALSE),
#                                         ifElse(identical(prepFun,cor),list(use = "pairwise.complete.obs"),list())),
#                          IsingLL = list(),
#                          huge = list(),
#                          adalasso = list()
#       )
#       
#       
#     }
#     
#     # estFun:
#     if (missing(estFun)){
#       estFun <- switch(default,
#                        EBICglasso = qgraph::EBICglasso,
#                        pcor = corpcor::cor2pcor,
#                        IsingFit = IsingFit::IsingFit,
#                        IsingLL = IsingSampler::EstimateIsing,
#                        huge = function(x)huge::huge.select(huge::huge(x,method = "glasso",verbose=FALSE), criterion = "ebic",verbose = FALSE),
#                        adalasso = parcor::adalasso.net
#       )
#     }
#     
#     # estArgs:
#     if (missing(estArgs)){
#       estArgs <- switch(default,
#                         EBICglasso = list(n = Np, returnAllResults = TRUE),
#                         IsingFit = list(plot = FALSE, progress = FALSE),
#                         pcor = list(),
#                         IsingLL = list(method = "ll"),
#                         huge = list(),
#                         adalasso = list()
#       )
#     }
#     
#     # graphFun:
#     if (missing(graphFun)){
#       graphFun <- switch(default,
#                          EBICglasso = function(x)x[['optnet']],
#                          IsingFit = function(x)x[['weiadj']],
#                          pcor = function(x)as.matrix(Matrix::forceSymmetric(x)),
#                          IsingLL = function(x)x[['graph']],
#                          huge = function(x)as.matrix(qgraph::wi2net(as.matrix(x$opt.icov))),
#                          adalasso = function(x)as.matrix(Matrix::forceSymmetric(x$pcor.adalasso))
#       )
#     }
#     
#     # graphArgs:
#     if (missing(graphArgs)){
#       graphArgs <- switch(default,
#                           EBICglasso = list(),
#                           IsingFit = list(),
#                           pcor = list(),
#                           IsingLL = list(),
#                           huge = list(),
#                           adalasso = list()
#       )
#     }
#     
#     # intFun:
#     if (missing(intFun)){
#       intFun <- switch(default,
#                        EBICglasso = null,
#                        IsingFit = function(x)x[['thresholds']],
#                        pcor = null,
#                        IsingLL = function(x) x[['thresholds']],
#                        huge = null,
#                        adalasso = null
#       )
#     }
#     
#     
#   }
#   
#   if (missing(prepFun)){
#     prepFun <- identity
#   }
#   
#   if (missing(prepArgs)){
#     prepArgs <- list()
#   }
#   
#   if (missing(graphFun)){
#     graphFun <- identity
#   }
#   
#   if (missing(graphArgs)){
#     graphArgs <- list()
#   }
#   
#   if (missing(intFun)){
#     intFun <- null
#   }
#   
#   if (missing(intArgs)){
#     intArgs <- list()
#   }
  
  # Compute args:
  args <- checkInput(
    default = default,
    prepFun = prepFun, # Fun to produce the correlation or covariance matrix
    prepArgs = prepArgs, # list with arguments for the correlation function
    estFun=estFun, # function that results in a network
    estArgs=estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
    graphFun=graphFun, # set to identity if missing
    graphArgs=graphArgs, # Set to null if missing
    intFun=intFun, # Set to null if missing
    intArgs=intArgs, # Set to null if missing
    sampleSize = Np
  )
  
  # Add data:
  args$data <- data

  # Compute input:
  input <- do.call(args$prepFun, c(list(data), args$prepArgs))
  
  # Compute network:
  res <- do.call(args$estFun, c(list(input),args$estArgs))
  
  
  sampleGraph <- do.call(args$graphFun,c(list(res), args$graphArgs))
  sampleResult <- list(
    graph = sampleGraph,
    intercepts = do.call(args$intFun,c(list(res), args$intArgs)),
    results = res,
    labels = labels,
    nNodes = ncol(data),
    nPerson = Np,
    input = args
#     input = list(
#       data = data,
#       default = default,
#       prepFun = prepFun, # Fun to produce the correlation or covariance matrix
#       prepArgs = prepArgs, # list with arguments for the correlation function
#       estFun = estFun, # function that results in a network
#       estArgs = estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
#       graphFun = graphFun, # set to identity if missing
#       graphArgs = graphArgs, # Set to null if missing
#       intFun = intFun, # Set to null if missing
#       intArgs = intArgs
#     )
  )
  class(sampleResult) <- c("bootnetResult", "list")
  
  # Plot?
#   if (plot){
#     plot(sampleResult,labels=labels,layout = "spring", ...)
#   }
  
  # Return network:
  return(sampleResult)
}