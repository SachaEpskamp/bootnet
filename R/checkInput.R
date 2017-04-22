
# Function that checks input and returns the functions:
checkInput <- function(
  default = c("none", "EBICglasso", "pcor","IsingFit","IsingSampler", "huge","adalasso","mgm","relimp"),
  fun, # Estimator function
  prepFun, # Fun to produce the correlation or covariance matrix
  prepArgs, # list with arguments for the correlation function
  estFun, # function that results in a network
  estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
  graphFun, # set to identity if missing
  graphArgs, # Set to null if missing
  intFun, # Set to null if missing
  intArgs, # Set to null if missing
  sampleSize,
  verbose=TRUE,
  construct = c("default","function","arguments"),
  .dots = list(),
  ... # Arguments to the estimator function
){
  if (default[[1]]=="glasso") default <- "EBICglasso"
  if (default[[1]]=="IsingSampler") default <- "IsingSampler"
  default <- match.arg(default)
  construct <- match.arg(construct)

  ### DEFAULT OPTIONS ###
 if (missing(fun)){
   fun <- NULL
 }
  

  
  # Stop if not compatible:
  dots <- c(.dots,list(...))
  
  # gather names:
  argNames <- character(0)
  
  if (!missing(prepFun)){
    argNames <- c(argNames,"prepFun")
  }
  if (!missing(prepArgs)){
    argNames <- c(argNames,"prepArgs")
  }
  if (!missing(estFun)){
    argNames <- c(argNames,"estFun")
  }
  if (!missing(estArgs)){
    argNames <- c(argNames,"estArgs")
  }
  if (!missing(graphFun)){
    argNames <- c(argNames,"graphFun")
  }
  if (!missing(graphArgs)){
    argNames <- c(argNames,"graphArgs")
  }
  if (!missing(intFun)){
    argNames <- c(argNames,"intFun")
  }
  if (!missing(intArgs)){
    argNames <- c(argNames,"intArgs")
  }

  # Not compatible if construct is used:
  if (length(dots) > 0 && construct == "arguments"){
    
    stop(paste0("Ambiguous argument specification. Old functonality is used (construct = 'arguments') in combination with new functionality arguments (implying construct = 'function'): ",
                paste0("'",names(dots),"'",collapse="; "),". These arguments are NOT compatible!"))
    
  }
  
  # relimp not compatable with old:
  if (construct == "arguments" & default == "relimp"){
    stop("default = 'relimp' not supported with old bootnet style (construct = 'arguments')")
    
  }
  
  if (length(argNames) > 0 && construct == "function"){
    
    stop(paste0("Ambiguous argument specification. New functonality is used (construct = 'function') in combination with old functionality arguments (implying construct = 'arguments'): ",
                paste0("'",argNames,"'",collapse="; "),". These arguments are NOT compatible!"))
    
  }
    
  # not compatible if both dots are used and arguments are used:
  if (length(argNames) > 0 & length(dots) > 0){

    stop(paste0("Ambiguous argument specification. Both old functionality arguments are used, compatible with construct = 'arguments': ",
                paste0("'",argNames,"'",collapse="; "),", as well as new functionality arguments are used, compatible with construct = 'function': ",
                paste0("'",names(dots),"'",collapse="; "),". These two types of arguments are NOT compatible!"))
    
  }
  
  
  # Check to construct via function or to construct via arguments:
  # if no default and no fun, use arguments:
  if (construct == "default"){
    construct <- "function"
    
    if (default == "none" && is.null(fun)){
      construct <- "arguments"
    }
    
    # If fun is missing, default is not none and one argument is not missing, use arguments (backward competability):
    if (default != "none" && is.null(fun) && (!missing(prepFun) | !missing(prepArgs) | !missing(estFun) | !missing(estArgs))){
      construct <- "arguments"
    }
  }
  
  # Check if arguments are not missing:
  if (default == "none" && construct == "arguments"){
    if (missing(prepFun) | missing(prepArgs) | missing(estFun) | missing(estArgs)){
      stop("If 'default' is not set and 'fun' is missing, 'prepFun', 'prepArgs', 'estFun' and 'estArgs' may not be missing.")
    }
  }
  
  ### Construct estimator function via function:
  if (construct == "function"){
    # Arguments:
    Args <- dots
    
    # Warn user that arguments are ignored:
    if (!missing(prepFun)){
      warning("'prepFun' argument is ignored as a function is used as arguments. To use 'prepFun', please set construct = 'arguments'")
    }
    if (!missing(prepArgs)){
      warning("'prepArgs' argument is ignored as a function is used as arguments. To use 'prepArgs', please set construct = 'arguments'")
    }
    if (!missing(estFun)){
      warning("'estFun' argument is ignored as a function is used as arguments. To use 'estFun', please set construct = 'arguments'")
    }
    if (!missing(estArgs)){
      warning("'estArgs' argument is ignored as a function is used as arguments. To use 'estArgs', please set construct = 'arguments'")
    }
    if (!missing(graphFun)){
      warning("'graphFun' argument is ignored as a function is used as arguments. To use 'graphFun', please set construct = 'arguments'")
    }
    if (!missing(graphArgs)){
      warning("'graphArgs' argument is ignored as a function is used as arguments. To use 'graphArgs', please set construct = 'arguments'")
    }
    if (!missing(intFun)){
      warning("'intFun' argument is ignored as a function is used as arguments. To use 'intFun', please set construct = 'arguments'")
    }
    if (!missing(intArgs)){
      warning("'intArgs' argument is ignored as a function is used as arguments. To use 'intArgs', please set construct = 'arguments'")
    }
    
    # per default:
    if (default == "none"){
      Function <- fun
    } else if (default == "EBICglasso"){
      Function <- bootnet_EBICglasso
    } else if (default == "IsingFit"){
      Function <- bootnet_IsingFit
    } else if (default == "IsingSampler"){
      Function <- bootnet_IsingSampler
    } else if (default == "pcor"){
      Function <- bootnet_pcor
    } else if (default == "adalasso"){
      Function <- bootnet_adalasso
    } else if (default == "huge"){
      Function <- bootnet_huge
    } else if (default == "mgm"){
      Function <- bootnet_mgm
    } else if (default == "relimp"){
      Function <- bootnet_relimp
    } else stop("Currently not supported.")
    
  } else {
    warning("Arguments (prepFun, estFun, etcetera) used to construct estimator. This functionality is deprecated and will no longer be supported in a future version of bootnet. Please consult the manual or contact the authors.")
    
    # Check dots, and warn user:
    if (length(dots) > 0){
      dotNames <- names(dots)
      warning(paste0("Arguments (prepFun, estFun, etcetera) used to construct estimator. As a result, the following arguments are ignored: ",paste0("'",dotNames,"'", collapse = ", "),". To use these arguments use construct = 'function' and supply a default set or set the 'fun' argument. In addition, do not use the 'prepFun', 'estFun', etcetera arguments."))
    }
    
    # Construct via arguments
    if (!(default == "none")){
      # prepFun:
      if (missing(prepFun)){
        prepFun <- switch(default,
                          EBICglasso = qgraph::cor_auto,
                          IsingFit = binarize,
                          IsingSampler = binarize,
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
                           IsingSampler = list(),
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
                         IsingSampler = IsingSampler::EstimateIsing,
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
                          IsingSampler = list(method = "ll"),
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
                           IsingSampler = function(x)x[['graph']],
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
                            IsingSampler = list(),
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
                         IsingSampler = function(x) x[['thresholds']],
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
    
    # Function:
    Function <- bootnet_argEstimator
    
    # List of arguents:
    Args <- list(
      prepFun = prepFun,
      prepArgs = prepArgs,
      estFun = estFun,
      estArgs = estArgs,
      graphFun = graphFun,
      graphArgs = graphArgs,
      intFun = intFun,
      intArgs = intArgs
    )
  }
  

  
 
  # Output:
  Output <- list(
    data = data,
    default = default,
    estimator = Function,
    arguments = Args
  )
  
  return(Output)
}