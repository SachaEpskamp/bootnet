# Null function:
null <- function(...) NULL

### IsingFit ###
# prep fun (make binary if needed):
binarize <- function(x, split = "median", na.rm=TRUE, removeNArows = TRUE, verbose = TRUE){
  x <- as.data.frame(x)
  
  if (all(unlist(x) %in% c(0,1))){
    return(x)
  } else {
    
    if (is.function(split)){
      splitName <- deparse(substitute(split))
    } else {
      splitName <- split
    }
    if (verbose){
      warning(paste("Splitting data by",splitName))
    }
    
    if (is.character(split) || is.function(split)){
      splits <- sapply(x,split,na.rm=na.rm)      
    } else {
      splits <- rep(split, length=ncol(x))
    }
    for (i in seq_len(ncol(x))){
      x[,i] <- 1 * (x[,i] < splits[i])
    }
    
    if (removeNArows){
      x <- x[apply(x,1,function(xx)all(!is.na(xx))),,drop=FALSE]
    }
    
    return(x)
  }
}


### ARGUMENT ESTIMATOR ###
# Construct the estimator:
bootnet_argEstimator <- function(data, prepFun, prepArgs, estFun, estArgs, graphFun, graphArgs, intFun, intArgs, verbose = TRUE){
  # prepArgs$verbose <- verbose
  if ("verbose" %in% names(formals(prepFun))){
    prepArgs$verbose <- verbose
  }
  
  # Compute input:
  input <- do.call(prepFun, c(list(data), prepArgs))
  
  # Compute network:
  res <- do.call(estFun, c(list(input),estArgs))
  
  # Extract network:
  sampleGraph <- do.call(graphFun,c(list(res), graphArgs))
  
  # Extract graph:
  intercepts <-  do.call(intFun,c(list(res), intArgs))
  
  # Return:
  return(list(graph=sampleGraph, intercepts = intercepts))
}

### EBIC GLASSO ESTIMATOR ###
bootnet_EBICglasso <- function(
  data, # Dataset used
  tuning = 0.5, # tuning parameter
  corMethod = c("cor_auto","cov","cor","npn"), # Correlation method
  missing = c("pairwise","listwise","stop"),
  sampleSize = c("maximum","minimim"), # Sample size when using missing = "pairwise"
  verbose = verbose,
  corArgs = list() # Extra arguments to the correlation function
){
  # Check arguments:
  corMethod <- match.arg(corMethod)
  missing <- match.arg(missing)
  sampleSize <- match.arg(sampleSize)
  
  # First test if data is a data frame:
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  # Obtain info from data:
  N <- ncol(data)
  Np <- nrow(data)
  
  # Check missing:
  if (missing == "stop"){
    if (any(is.na(data))){
      stop("Missing data detected and missing = 'stop'")
    }
  }
  
  # Correlate data:
  # npn:
  if (corMethod == "npn"){
    data <- huge::huge.npn(data)
    corMethod <- "cor"
  }
  
  # cor_auto:
  if (corMethod == "cor_auto"){
    args <- list(data=data,missing=missing,verbose=verbose)
    if (length(corArgs) > 0){
      for (i in seq_along(corArgs)){
        args[[names(corArgs)[[i]]]] <- corArgs[[i]]
      }
    }
    corMat <- do.call(qgraph::cor_auto,args)
  } else if (corMethod%in%c("cor","cov")){
    # Normal correlations
    
    use <- switch(missing,
                  pairwise = "pairwise.complete.obs",
                  listwise = "complete.obs")
    
    args <- list(x=data,use=use)
    if (length(corArgs) > 0){
      for (i in seq_along(corArgs)){
        args[[names(corArgs)[[i]]]] <- corArgs[[i]]
      }
    }
    
    corMat <- do.call(corMethod,args)
  } else stop ("Correlation method is not supported.")
  
  # Sample size:
  if (missing == "listwise"){
    sampleSize <- nrow(na.omit(data))
  } else{
    if (sampleSize == "maximum"){
      sampleSize <- sum(apply(data,1,function(x)!all(is.na(x))))
    } else {
      sampleSize <- sum(apply(data,1,function(x)!any(is.na(x))))
    }
  } 
  
  # Estimate network:
  Results <- qgraph::EBICglasso(corMat,
                                n =  sampleSize, 
                                gamma = tuning,
                                returnAllResults = TRUE)
  
  # Return:
  return(list(graph=Results$optnet,results=Results))
}


### PCOR ESTIMATOR ###
bootnet_pcor <- function(
  data, # Dataset used
  corMethod = c("cor_auto","cov","cor","npn"), # Correlation method
  missing = c("pairwise","listwise","stop"),
  verbose = verbose,
  corArgs = list() # Extra arguments to the correlation function
){
  # Check arguments:
  corMethod <- match.arg(corMethod)
  missing <- match.arg(missing)
  
  # First test if data is a data frame:
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  # Obtain info from data:
  N <- ncol(data)
  Np <- nrow(data)
  
  # Check missing:
  if (missing == "stop"){
    if (any(is.na(data))){
      stop("Missing data detected and missing = 'stop'")
    }
  }
  
  # Correlate data:
  # npn:
  if (corMethod == "npn"){
    data <- huge::huge.npn(data)
    corMethod <- "cor"
  }
  
  # cor_auto:
  if (corMethod == "cor_auto"){
    args <- list(data=data,missing=missing,verbose=verbose)
    if (length(corArgs) > 0){
      for (i in seq_along(corArgs)){
        args[[names(corArgs)[[i]]]] <- corArgs[[i]]
      }
    }
    corMat <- do.call(qgraph::cor_auto,args)
  } else if (corMethod%in%c("cor","cov")){
    # Normal correlations
    
    use <- switch(missing,
                  pairwise = "pairwise.complete.obs",
                  listwise = "complete.obs")
    
    args <- list(x=data,use=use)
    if (length(corArgs) > 0){
      for (i in seq_along(corArgs)){
        args[[names(corArgs)[[i]]]] <- corArgs[[i]]
      }
    }
    
    corMat <- do.call(corMethod,args)
  } else stop ("Correlation method is not supported.")
  
  
  # Estimate network:
  Results <- getWmat(qgraph::qgraph(corMat,graph = "pcor",DoNotPlot = TRUE))
  
  # Return:
  return(list(graph=Results,results=Results))
}


### ISINGFIT ESTIMATOR ###
bootnet_IsingFit <- function(
  data, # Dataset used
  tuning = 0.25, # tuning parameter
  missing = c("listwise","stop"),
  verbose = verbose,
  rule = c("AND","OR"),
  split = "median"
){
  # Check arguments:
  missing <- match.arg(missing)
  rule <- match.arg(rule)
  
  # First test if data is a data frame:
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  
  # Check missing:
  if (missing == "stop"){
    if (any(is.na(data))){
      stop("Missing data detected and missing = 'stop'")
    }
  } else {
    # listwise:
    data <- na.omit(data)
  }
  
  
  # Obtain info from data:
  N <- ncol(data)
  Np <- nrow(data)
  
  # Binarize:
  if (!all(data %in% c(0,1))){
    data <- bootnet::binarize(data, split = split, verbose = verbose)
  }
  
  # Estimate network:
  Results <- IsingFit::IsingFit(data, AND = rule == "AND", gamma = tuning,progressbar = verbose,plot = FALSE)
  
  # Return:
  return(list(graph = Results$weiadj, intercepts = Results$thresholds,
              results = Results))
}

### ISINGSAMPLER ESTIMATOR ###
bootnet_IsingSampler <- function(
  data, # Dataset used
  missing = c("listwise","stop"),
  verbose = verbose,
  split = "median",
  method = c("default","ll","pl","uni","bi")
){
  # Check arguments:
  missing <- match.arg(missing)
  method <- match.arg(method)
  
  # First test if data is a data frame:
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  if (method == "default"){
    if (ncol(data) > 20){
      method <- "uni"
      if (verbose){
        message("'method' set to 'uni'")
      }
    } else {
      method <- "ll"
      if (verbose){
        message("'method' set to 'll'")
      }
    }
  }
  
  # Check missing:
  if (missing == "stop"){
    if (any(is.na(data))){
      stop("Missing data detected and missing = 'stop'")
    }
  } else {
    # listwise:
    data <- na.omit(data)
  }
  
  
  # Obtain info from data:
  N <- ncol(data)
  Np <- nrow(data)
  
  # Binarize:
  if (!all(data %in% c(0,1))){
    data <- bootnet::binarize(data, split = split, verbose = verbose)
  }
  
  # Estimate network:
  Results <- IsingSampler::EstimateIsing(as.matrix(data), method = method)
  
  # Return:
  return(list(graph = Results$graph, intercepts = Results$thresholds,
              results = Results))
}


### PCOR ESTIMATOR ###
bootnet_adalasso <- function(
  data, # Dataset used
  missing = c("listwise","stop"),
  verbose = verbose,
  nFolds = 10 # Number of folds
){
  # Check arguments:
  missing <- match.arg(missing)
  
  # First test if data is a data frame:
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  # Obtain info from data:
  N <- ncol(data)
  Np <- nrow(data)
  
  
  # Check missing:
  if (missing == "stop"){
    if (any(is.na(data))){
      stop("Missing data detected and missing = 'stop'")
    }
  } else {
    # listwise:
    data <- na.omit(data)
  }
  
  # Estimate network:
  Results <- parcor::adalasso.net(data, k = nFolds)
  
  # Return:
  return(list(graph=as.matrix(Matrix::forceSymmetric(Results$pcor.adalasso)),results=Results))
}


### HUGE ESTIMATOR ###
bootnet_huge <- function(
  data, # Dataset used
  tuning = 0.5,
  missing = c("listwise","stop"),
  verbose = verbose,
  npn = TRUE, # Compute nonparanormal?
  criterion = c("ebic","ric","stars")
  # method = c("glasso","mb","ct")
){
  # Check arguments:
  missing <- match.arg(missing)
  criterion <- match.arg(criterion)
  # method <- match.arg(method)
  
  # First test if data is a data frame:
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  # Obtain info from data:
  N <- ncol(data)
  Np <- nrow(data)
  
  
  # Check missing:
  if (missing == "stop"){
    if (any(is.na(data))){
      stop("Missing data detected and missing = 'stop'")
    }
  } else {
    # listwise:
    data <- na.omit(data)
  }
  
  # Nonparanormal:
  if (npn){
    data <- huge::huge.npn(na.omit(as.matrix(data)),verbose = verbose)
  }
  
  # Estimate network:
  Results <- huge::huge.select(huge::huge(data,method = "glasso",verbose=verbose), criterion = criterion,verbose = verbose,ebic.gamma = tuning)
  
  # Return:
  return(list(
    graph=as.matrix(qgraph::wi2net(as.matrix(Results$opt.icov))),
    results=Results))
}

### MGM ESTIMATOR ###
bootnet_mgm <- function(
  data, # Dataset used
  type,
  lev,
  tuning = 0.5,
  missing = c("listwise","stop"),
  verbose = verbose,
  criterion = c("EBIC","CV"),
  nFolds = 10,
  degree = 2,
  rule = c("AND","OR")
  # method = c("glasso","mb","ct")
){
  # Check arguments:
  missing <- match.arg(missing)
  criterion <- match.arg(criterion)
  # method <- match.arg(method)
  
  # First test if data is a data frame:
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  # Obtain info from data:
  N <- ncol(data)
  Np <- nrow(data)
  
  
  # Check missing:
  if (missing == "stop"){
    if (any(is.na(data))){
      stop("Missing data detected and missing = 'stop'")
    }
  } else {
    # listwise:
    data <- na.omit(data)
  }
  
  # Set type automatically:
  if (missing(type)){
    if (verbose){
      message("'type' argument not assigned. Setting type to 'c' for all binary variables and 'g' for all other variables.")     
    }
    type <- ifelse(apply(data,2,function(x)all(x%in%c(0,1))),"c","g")
  }
  
  # Set lev automatically:
  if (missing(lev)){
    if (verbose){
      message("'lev' argument not assigned. Setting lev to 1 for all Gaussian/Poisson variables and number of unique values for all categorical variables")  
    }
    
    lev <- ifelse(type == "c", apply(data,2,function(x)length(unique(x))),1)
  }
  
  # Estimate:
  log <- capture.output(Results <- mgmfit(data,type=type,lev=lev,lambda.sel = criterion,folds = nFolds,gam = tuning,d = degree,pbar = verbose))
  
  # Warn for unsigned:
  if (any(Results$signs==0,na.rm = TRUE)){
    warning("Bootnet does not support unsigned edges and treats these as positive edges.")
  }
  Results$signs[is.na(Results$signs)] <- 0
  
  # Graph:
  Graph <- Results$wadj
  Graph <- ifelse(Results$signs==-1,-Graph,Graph)
  
  # Return:
  return(list(
    graph=Graph,
    results=Results))
}