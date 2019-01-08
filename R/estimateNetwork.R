# This function takes data as input and produced a network. It is used inside bootnet:
estimateNetwork <- function(
  data,
  default = c("none", "EBICglasso", "pcor","IsingFit","IsingSampler", "huge","adalasso","mgm","relimp", "cor","TMFG",
              "ggmModSelect", "LoGo","graphicalVAR", "piecewiseIsing","SVAR_lavaan"),
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
  datatype,
  checkNumeric = FALSE,
  # plot = TRUE, # Plot the network?
  ..., # Arguments to the 'fun' function
  .input, # Skips most of first steps if supplied
  memorysaver = FALSE # If set to FALSE data, estimator and results are not stored.
){
  # Borsboom easter egg:
  if (default[1] == "Borsboom") return(42)
  
  if (default[[1]]=="glasso") default <- "EBICglasso"
  default <- match.arg(default)
  
  # datatype test:
  if (missing(datatype)){
    if (is(data,"tsData")){
      datatype <- "graphicalVAR"
    } else {
      datatype <- "normal"
    }
      
  }

  if (!datatype%in% c("normal","graphicalVAR")){
    stop("Only datatypes 'normal' and 'graphicalVAR' currently supported.")
  }
#   
  # If NAs and default can't handle, stop:
  # if (any(is.na(data)) && default %in% c("huge","adalasso")){
  #   stop(paste0("Missing data not supported for default set '",default,"'. Try using na.omit(data)."))
  # }

  # First test if data is a data frame:
  if (datatype == "normal" && !(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (datatype == "normal" && is.matrix(data)){
    data <- as.data.frame(data)
  }
 
  if (missing(directed)){
    if (default == "graphicalVAR"){
      directed <- list(contemporaneous = FALSE, temporal = TRUE)
    } else  if (default == "SVAR_lavaan"){
      directed <- list(contemporaneous = TRUE, temporal = TRUE)
    } else if (!default %in% c("relimp","DAG")){
      directed <- FALSE 
    } else {
      directed <- TRUE
    }
  }
  
  if (datatype == "normal"){
    N <- ncol(data)
    Np <- nrow(data)   
    if (missing(labels)){
      labels <- colnames(data)
    }
    
    if (checkNumeric){
      # Check and remove any variable that is not ordered, integer or numeric:
      goodColumns <- sapply(data, function(x) is.numeric(x) | is.ordered(x) | is.integer(x))
      
      if (!all(goodColumns)){
        if (verbose){
          warning(paste0("Removing non-numeric columns: ",paste(which(!goodColumns),collapse="; ")))
        }
        data <- data[,goodColumns,drop=FALSE]
      }
    }

    
  } else if (datatype == "graphicalVAR"){
   N <- length(data$vars)
   Np <- nrow(data$data_c)
   if (missing(labels)){
     labels <- data$vars
   }
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

  if (!is.list(Result)){
    sampleGraph <- Result
    intercepts <- NULL
    output <- Result
    nNode <- ncol(Result)
  } else if (is.list(Result$graph)){
    sampleGraph <- Result$graph
    intercepts <- Result$intercepts
    output <- Result$results
    nNode <- ncol(Result$graph[[1]])
  } else {
    sampleGraph <- Result$graph
    intercepts <- Result$intercepts
    output <- Result$results
    nNode <- ncol(Result$graph)
  }

  
  if (!is.matrix(sampleGraph)){
    if (is.list(sampleGraph)){
      if (!is.matrix(sampleGraph[[1]])){
        stop("Estimated result is not a list of matrices encoding networks.")   
      }
    } else {
      stop("Estimated result is not a matrix encoding a network.")
    }
  }
  
  # Special data?
  if (!is.list(Result) || is.null(Result$specialData)){
    outdata <- data
    datatype <- "normal"
  } else {
    outdata <- Result$specialData$data
    datatype <- Result$specialData$type
  }
  

  sampleResult <- list(
    graph = sampleGraph,
    intercepts = intercepts,
    results = output,
    labels = labels,
    nNode = nNode,
    nPerson = Np,
    estimator = .input$estimator,
    arguments = .input$arguments,
    data = outdata,
    datatype = datatype,
    default = default,
    weighted = weighted,
    signed = signed,
    directed=directed,
    .input = .input,
    thresholded = FALSE
  )
  class(sampleResult) <- c("bootnetResult", "list")
  
  if (default == "graphicalVAR"){
       sampleResult$labels <- output$data$vars
  }
  if (default == "SVAR_lavaan"){
    sampleResult$labels  <- outdata$vars
  }
  
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