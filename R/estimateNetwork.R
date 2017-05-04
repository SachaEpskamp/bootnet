# This function takes data as input and produced a network. It is used inside bootnet:
estimateNetwork <- function(
  data,
  default = c("none", "EBICglasso", "pcor","IsingFit","IsingSampler", "huge","adalasso","mgm","relimp"),
  fun, # A function that takes data and returns a network or list entitled "graph" and "thresholds". optional.
  prepFun, # Fun to produce the correlation or covariance matrix
  prepArgs, # list with arguments for the correlation function
  estFun, # function that results in a network
  estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
  graphFun, # set to identity if missing
  graphArgs, # Set to null if missing
  intFun, # Set to null if missing
  intArgs, # Set to null if missing
  labels, # if missing taken from colnames
  verbose = TRUE, # Dummy used in cor_auto and in the future in other functions. Set to FALSE in bootnet
  construct = c("default","function","arguments"),
  .dots = list(),
  weighted = TRUE,
  signed = TRUE,
  directed,
  # plot = TRUE, # Plot the network?
  ..., # Arguments to the 'fun' function
  .input, # Skips most of first steps if supplied
  memorysaver = FALSE # If set to FALSE data, estimator and results are not stored.
){
  if (default[[1]]=="glasso") default <- "EBICglasso"
  default <- match.arg(default)
#   
  # If NAs and default can't handle, stop:
  # if (any(is.na(data)) && default %in% c("huge","adalasso")){
  #   stop(paste0("Missing data not supported for default set '",default,"'. Try using na.omit(data)."))
  # }

  # First test if data is a data frame:
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
 
  if (missing(directed)){
    if (!default %in% c("graphicalVAR","relimp","DAG")){
      directed <- FALSE 
    } else {
      directed <- TRUE
    }
  }
  
  N <- ncol(data)
  Np <- nrow(data)
  
  
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
  
  
  # Compute estimator:
  if (missing(.input)){
    .input <- checkInput(
      default = default,
      fun = fun,
      prepFun = prepFun, # Fun to produce the correlation or covariance matrix
      prepArgs = prepArgs, # list with arguments for the correlation function
      estFun=estFun, # function that results in a network
      estArgs=estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
      graphFun=graphFun, # set to identity if missing
      graphArgs=graphArgs, # Set to null if missing
      intFun=intFun, # Set to null if missing
      intArgs=intArgs, # Set to null if missing
      sampleSize = Np,
      construct=construct,
      verbose=verbose,
      .dots=.dots,
      ...
    )
  }

  
  # Add verbose:
  # Every estimator must have argument verbose:
  if ("verbose" %in% names(formals(.input$estimator))){
    .input$arguments$verbose <- verbose
  }

  # Compute network:
  Result <- do.call(.input$estimator, c(list(data),.input$arguments))
  
  if (is.list(Result)){
    sampleGraph <- Result$graph
    intercepts <- Result$intercepts
    output <- Result$results
  } else {
    sampleGraph <- Result
    intercepts <- NULL
    output <- NULL
  }
  
  if (!is.matrix(sampleGraph)){
    stop("Estimated result is not a matrix encoding a network.")
  }

  sampleResult <- list(
    graph = sampleGraph,
    intercepts = intercepts,
    results = output,
    labels = labels,
    nNodes = ncol(data),
    nPerson = Np,
    estimator = .input$estimator,
    arguments = .input$arguments,
    data = data,
    default = default,
    weighted = weighted,
    signed = signed,
    directed=directed,
    .input = .input
  )
  class(sampleResult) <- c("bootnetResult", "list")
  
  # Memory save:
  if(memorysaver)
  {
    sampleResult$results <- NA
    sampleResult$estimator <- NA
    sampleResult$data <- NA
    sampleResult$.input <- NA
  }
  
  # Plot?
#   if (plot){
#     plot(sampleResult,labels=labels,layout = "spring", ...)
#   }
  
  # Return network:
  return(sampleResult)
}