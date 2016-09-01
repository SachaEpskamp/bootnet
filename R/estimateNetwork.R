# This function takes data as input and produced a network. It is used inside bootnet:
estimateNetwork <- function(
  data,
  default = c("none", "EBICglasso", "pcor","IsingFit","IsingLL", "huge","adalasso"),
  prepFun, # Fun to produce the correlation or covariance matrix
  prepArgs, # list with arguments for the correlation function
  estFun, # function that results in a network
  estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
  graphFun, # set to identity if missing
  graphArgs, # Set to null if missing
  intFun, # Set to null if missing
  intArgs, # Set to null if missing
  labels, # if missing taken from colnames
  verbose = TRUE # Dummy used in cor_auto and in the future in other functions. Set to FALSE in bootnet
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
    sampleSize = Np,
    verbose=verbose
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