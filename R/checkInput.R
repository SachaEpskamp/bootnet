
# Function that checks input and returns the functions:
checkInput <- function(
  default = c("none", "EBICglasso", "pcor","IsingFit","IsingLL", "huge","adalasso"),
  prepFun, # Fun to produce the correlation or covariance matrix
  prepArgs, # list with arguments for the correlation function
  estFun, # function that results in a network
  estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
  graphFun, # set to identity if missing
  graphArgs, # Set to null if missing
  intFun, # Set to null if missing
  intArgs, # Set to null if missing
  sampleSize,
  verbose=TRUE
){
  if (default[[1]]=="glasso") default <- "EBICglasso"
  default <- match.arg(default)
  
  ### DEFAULT OPTIONS ###
  if ((default == "none")){
    if (missing(prepFun) | missing(prepArgs) | missing(estFun) | missing(estArgs)){
      stop("If 'default' is not set, 'prepFun', 'prepArgs', 'estFun' and 'estArgs' may not be missing.")
    }
  }
  
  if (!(default == "none")){
    # prepFun:
    if (missing(prepFun)){
      prepFun <- switch(default,
                        EBICglasso = qgraph::cor_auto,
                        IsingFit = binarize,
                        IsingLL = binarize,
                        pcor = qgraph::cor_auto,
                        huge = function(x)huge::huge.npn(na.omit(as.matrix(x)),verbose = FALSE),
                        adalasso = identity
      )
      #       prepFun <- switch(default,
      #                         EBICglasso = cor,
      #                         IsingFit = binarize,
      #                         pcor = cor
      #       )      
    }
    
    # prepArgs:
    #     qgraphVersion <- packageDescription("qgraph")$Version
    #     qgraphVersion <- as.numeric(strsplit(qgraphVersion,split="\\.|\\-")[[1]])
    #     if (length(qgraphVersion)==1) qgraphVersion <- c(qgraphVersion,0)
    #     if (length(qgraphVersion)==2) qgraphVersion <- c(qgraphVersion,0)
    #     goodVersion <- 
    #       (qgraphVersion[[1]] >= 1 & qgraphVersion[[2]] >= 3 & qgraphVersion[[3]] >= 1) | 
    #       (qgraphVersion[[1]] >= 1 & qgraphVersion[[2]] > 3) | 
    #       qgraphVersion[[1]] > 1
    
    if (missing(prepArgs)){
      prepArgs <- switch(default,
                         EBICglasso = ifElse(identical(prepFun,qgraph::cor_auto),list(verbose=verbose),
                                             ifElse(identical(prepFun,cor),list(use = "pairwise.complete.obs"),list())),
                         IsingFit = list(),
                         pcor =  ifElse(identical(prepFun,qgraph::cor_auto),list(verbose=verbose),
                                        ifElse(identical(prepFun,cor),list(use = "pairwise.complete.obs"),list())),
                         IsingLL = list(),
                         huge = list(),
                         adalasso = list()
      )


    }

    
    # estFun:
    if (missing(estFun)){
      estFun <- switch(default,
                       EBICglasso = qgraph::EBICglasso,
                       pcor = corpcor::cor2pcor,
                       IsingFit = IsingFit::IsingFit,
                       IsingLL = IsingSampler::EstimateIsing,
                       huge = function(x)huge::huge.select(huge::huge(x,method = "glasso",verbose=FALSE), criterion = "ebic",verbose = FALSE),
                       adalasso = parcor::adalasso.net
      )
    }
    
    # estArgs:
    if (missing(estArgs)){
      estArgs <- switch(default,
                        EBICglasso = list(n = sampleSize, returnAllResults = TRUE),
                        IsingFit = list(plot = FALSE, progress = FALSE),
                        pcor = list(),
                        IsingLL = list(method = "ll"),
                        huge = list(),
                        adalasso = list()
      )
    }
    
    # graphFun:
    if (missing(graphFun)){
      graphFun <- switch(default,
                         EBICglasso = function(x)x[['optnet']],
                         IsingFit = function(x)x[['weiadj']],
                         pcor = function(x)as.matrix(Matrix::forceSymmetric(x)),
                         IsingLL = function(x)x[['graph']],
                         huge = function(x)as.matrix(qgraph::wi2net(as.matrix(x$opt.icov))),
                         adalasso = function(x)as.matrix(Matrix::forceSymmetric(x$pcor.adalasso))
      )
    }
    
    # graphArgs:
    if (missing(graphArgs)){
      graphArgs <- switch(default,
                          EBICglasso = list(),
                          IsingFit = list(),
                          pcor = list(),
                          IsingLL = list(),
                          huge = list(),
                          adalasso = list()
      )
    }
    
    # intFun:
    if (missing(intFun)){
      intFun <- switch(default,
                       EBICglasso = null,
                       IsingFit = function(x)x[['thresholds']],
                       pcor = null,
                       IsingLL = function(x) x[['thresholds']],
                       huge = null,
                       adalasso = null
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
  
  Args <- list(
    data = data,
    default = default,
    prepFun = prepFun, # Fun to produce the correlation or covariance matrix
    prepArgs = prepArgs, # list with arguments for the correlation function
    estFun = estFun, # function that results in a network
    estArgs = estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
    graphFun = graphFun, # set to identity if missing
    graphArgs = graphArgs, # Set to null if missing
    intFun = intFun, # Set to null if missing
    intArgs = intArgs
  )
  
  return(Args)
}