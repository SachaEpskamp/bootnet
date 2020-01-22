ifElse <- function(statement,true,false){
  if (statement){
    return(true)
  } else {
    return(false)
  }
}

noDiag <- function(x){
  diag(x) <- 0
  return(x)
}

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
  default = c("none", "EBICglasso", "ggmModSelect", "pcor","IsingFit","IsingSampler", "huge","adalasso","mgm","relimp","cor","TMFG", "ggmModSelect", "LoGo","SVAR_lavaan"), # Default method to use. EBICglasso, IsingFit, concentration, some more....
  type = c("nonparametric","parametric","node","person","jackknife","case"), # Bootstrap method to use
  nCores = 1,
  statistics = c("edge","strength","outStrength","inStrength"),
  model = c("detect","GGM","Ising","graphicalVAR"), # Models to use for bootstrap method = parametric. Detect will use the default set and estimation function.
  fun,
  # prepFun, # Fun to produce the correlation or covariance matrix
  # prepArgs, # list with arguments for the correlation function
  # estFun, # function that results in a network
  # estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
  # graphFun, # set to identity if missing
  # graphArgs, # Set to null if missing
  # intFun, # Set to null if missing
  # intArgs, # Set to null if missing
  verbose = TRUE, # messages on what is being done?
  # construct = c("default","function","arguments"),
  labels, # if missing taken from colnames
  alpha = 1, # centrality alpha
  caseMin = 0.05, # minimum proportion to DROP
  caseMax = 0.75, # Maximum proportion to DROP
  caseN = 10,
  subNodes, # if type = "node", defaults to 2:(p-1)
  subCases,
  computeCentrality = TRUE,
  propBoot = 1, # M out of N
  # subsampleSize,
  replacement = TRUE,
  graph, # for parametric bootstrap
  sampleSize, # for parametric bootstrap
  intercepts, # for parametric bootstrap
  weighted,
  signed,
  directed,
  includeDiagonal = FALSE,
  communities=NULL,
  useCommunities="all",
  library = .libPaths(),
  memorysaver = TRUE,
  # datatype = c("normal","graphicalVAR"), # Extracted from object or given
  ... # Other arguments
  # edgeResample = FALSE # If true, only resample edges from original estimate
  # scaleAdjust = FALSE
){
  construct <- "function"
  if (default[[1]]=="glasso") default <- "EBICglasso"
  default <- match.arg(default)
  
  # Check default:
  if (default == "graphicalVAR" && !is(data,"bootnetResult")){
    stop("default = 'graphicalVAR' only supported for output of estimateNetwork()")
  }
  
  # Check if statistics is all:
  if (any(statistics=="all")){
    statistics <- c("intercept","edge","length","distance","closeness","betweenness","strength","expectedInfluence",
                    "outStrength","outExpectedInfluence","inStrength","inExpectedInfluence","rspbc","hybrid",
                    "bridgeStrength", "bridgeCloseness", "bridgeBetweenness",
                    "bridgeExpectedInfluence")
  } else {
    message(paste("Note: bootnet will store only the following statistics: ",paste0(statistics, collapse=", ")))
  }
  
  
  type <- match.arg(type)
  # datatype <- match.arg(datatype)
  if (type == "case") type <- "person"
  model <- match.arg(model)
  
  # If data is bootnetResult, extract:
  
  # If data is missing, checks for parametric bootstrap:
  if (missing(data)){
    if (type != "parametric"){
      warning("'data' can only be missing if type = 'parametric'. Setting type = 'parametric' and performing parametric bootstrap instead.")
      type <- 'parametric'
    }
    
    if (missing(graph)){
      stop("'graph' may not be missing in parametric bootstrap when 'data' is missing.")
    }
    if (missing(sampleSize)){
      stop("'sampleSize' may not be missing in parametric bootstrap when 'data' is missing.")
    }
    
    
    N <- ncol(graph)
    Np <-  sampleSize
    datatype <- "normal"
    
    if (missing(intercepts)){
      intercepts <- rep(0, N)
    }
    
    if (!missing(data)){
      warning("'data' is ignored when using manual parametric bootstrap.")
      data <- NULL
    }
    manual <- TRUE
    dots <- list(...)
  } else {
    
    manual <- FALSE
    
    if (is(data,"bootnetResult")){
      
      # Check if thresholded:
      if (isTRUE(data$thresholded)){
        stop("Network has already been thresholded using bootstraps.")
      }
      
      # Check if bootInclude:
      if (isTRUE(data$bootInclude)){
        stop("Network is based on bootstrap include probabilities.")
      }
      
      default <- data$default
      inputCheck <- data$.input
      datatype <- data$datatype
      if (missing(labels)){
        labels <- data$labels
      }
      # prepFun <- data$input$prepFun
      # prepArgs <- data$input$prepArgs
      # estFun <- data$input$estFun
      # estArgs <- data$input$estArgs
      # graphFun <- data$input$graphFun
      # graphArgs <- data$input$graphArgs
      # intFun <- data$input$intFun
      # intArgs <- data$input$intArgs
      fun <- data$estimator
      dots <- data$arguments
      if (missing(weighted)){
        weighted <- data$weighted
      }
      if (missing(signed)){
        signed <- data$signed
      }
      if (missing(directed)){
        directed <- data$directed
      }
      
      N <- data$nNode
      Np <- data$nPerson
      data <- data$data
      
      
    } else {
      
      datatype <- "normal"
      dots <- list(...)
      N <- ncol(data)
      Np <- nrow(data)
      if (missing(fun)){
        fun <- NULL        
      }
      
      # Check and remove any variable that is not ordered, integer or numeric:
      if (!manual){
        goodColumns <- sapply(data, function(x) is.numeric(x) | is.ordered(x) | is.integer(x))
        
        if (!all(goodColumns)){
          if (verbose){
            warning(paste0("Removing non-numeric columns: ",paste(which(!goodColumns),collapse="; ")))
          }
          data <- data[,goodColumns,drop=FALSE]
        }
      }
      
      # 
      # inputCheck <- checkInput(
      #   default = default,
      #   fun = fun,
      #   prepFun = prepFun, # Fun to produce the correlation or covariance matrix
      #   prepArgs = prepArgs, # list with arguments for the correlation function
      #   estFun=estFun, # function that results in a network
      #   estArgs=estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
      #   graphFun=graphFun, # set to identity if missing
      #   graphArgs=graphArgs, # Set to null if missing
      #   intFun=intFun, # Set to null if missing
      #   intArgs=intArgs, # Set to null if missing
      #   sampleSize = Np,
      #   construct = construct,
      #   .dots = dots
      # )
      # 
      
    }
    
    
    
    
  }
  
  # subNodes and subCases:
  if (missing(subNodes)){
    if (datatype == "normal"){
      subNodes <- 2:(N-1)
    } else if (datatype == "graphicalVAR"){
      subNodes <- 2:(length(vars)-1)
    }
  }
  
  if (missing(subCases)){
    if (datatype == "normal"){
      subCases <- round((1-seq(caseMin,caseMax,length=caseN)) * Np)
    } else if (datatype == "graphicalVAR"){
      subCases <- round((1-seq(caseMin,caseMax,length=caseN)) * nrow(data$data_c))
    }
  }
  
  
  inputCheck <- checkInput(
    default = default,
    fun = fun,
    # prepFun = prepFun, # Fun to produce the correlation or covariance matrix
    # prepArgs = prepArgs, # list with arguments for the correlation function
    # estFun=estFun, # function that results in a network
    # estArgs=estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
    # graphFun=graphFun, # set to identity if missing
    # graphArgs=graphArgs, # Set to null if missing
    # intFun=intFun, # Set to null if missing
    # intArgs=intArgs, # Set to null if missing
    # sampleSize = Np,
    # construct = construct,
    .dots = dots
  )
  
  # Weighted and signed defaults
  if (missing(weighted)){
    weighted <- TRUE
  }
  if (missing(signed)){
    signed <- TRUE
  }
  if (missing(directed)){
    if (!default %in% c("graphicalVAR","relimp","DAG")) directed <- FALSE 
  }
  
  
  if (type == "jackknife"){
    message("Jacknife overwrites nBoot to sample size")
    nBoots <- Np
  }
  
  
  if (type == "node" & N < 3){
    stop("Node-wise bootstrapping requires at least three nodes.")
  }
  
  # First test if data is a data frame:
  if (datatype == "normal" && !manual && !(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If matrix coerce to data frame:
  if (!manual && is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  if (missing(labels)){
    if (manual){
      labels <- colnames(graph)
      if (is.null(labels)){
        labels <- seq_len(ncol(graph))
      }
    } else {
      labels <- colnames(data)    
      if (is.null(labels)){
        labels <- seq_len(ncol(data))
      }
    }
    
    
  }
  
  
  
  ## For parametric bootstrap, detect model
  if (type == "parametric" & model == "detect"){
    if (manual){
      stop("'model' must be set in parametric bootstrap without 'data'.")
    }
    
    if (default == "graphicalVAR") {
      model <- "graphicalVAR"  
    } else if (default != "none"){
      model <- ifelse(grepl("ising",default,ignore.case=TRUE),"Ising","GGM")
    } else {
      stop("'none' default set not supported for graphicalVAR data.")
      # model <- ifelse(any(grepl("ising",deparse(estFun),ignore.case=TRUE)),"Ising","GGM")
    }
    message(paste0("model set to '",model,"'"))
  }
  
  # Estimate sample result:
  # Check the input:
  
  
  if (!manual)
  {
    if (verbose){
      message("Estimating sample network...")
    }
    
    sampleResult <- estimateNetwork(data, 
                                    default = default,
                                    fun = inputCheck$estimator,
                                    .dots = inputCheck$arguments,
                                    labels = labels,
                                    verbose = verbose,
                                    weighted = weighted,
                                    signed = signed,
                                    .input = inputCheck,
                                    datatype = datatype)
    
    
  } else {
    
    sampleResult <- list(
      graph = graph,
      intercepts = intercepts,
      labels = labels,
      nNodes = N,
      nPerson = Np,
      estimator = inputCheck$estimator,
      arguments = inputCheck$arguments,
      default = default,
      weighted = weighted,
      signed = signed
    )
    class(sampleResult) <- c("bootnetResult", "list")
    
  }
  
  # Extract arguments:
  # default <- sampleResult$input$default
  # prepFun <- sampleResult$input$prepFun
  # prepArgs <- sampleResult$input$prepArgs
  # estFun <- sampleResult$input$estFun
  # estArgs <- sampleResult$input$estArgs
  # graphFun <- sampleResult$input$graphFun
  # graphArgs <- sampleResult$input$graphArgs
  # intFun <- sampleResult$input$intFun
  # intArgs <- sampleResult$input$intArgs
  
  
  # if (!isSymmetric(as.matrix(sampleResult[['graph']]))){
  #   stop("bootnet does not support directed graphs")
  # }
  
  
  #   ### Observation-wise bootstrapping!
  #   if (type == "observation"){
  # Bootstrap results:
  if (nCores == 1){
    bootResults <- vector("list", nBoots)
    
    if (verbose){
      message("Bootstrapping...")
      pb <- txtProgressBar(0,nBoots,style = 3)
    }
    
    for (b in seq_len(nBoots)){
      
      tryLimit <- 10
      tryCount <- 0
      repeat{
        
        if (! type %in% c("node","person")){
          nNode <- N
          inSample <- seq_len(N)
          
          if (type == "jackknife"){
            if (datatype == "normal"){
              bootData <- data[-b,,drop=FALSE]                  
            } else {
              bootData <- data
              bootData$data_c <- bootData$data_c[-b,,drop=FALSE]
              bootData$data_l <- bootData$data_l[-b,,drop=FALSE]    
            }
            
            nPerson <- Np - 1
          } else if (type == "parametric"){
            nPerson <- Np
            if (model == "Ising"){
              bootData <- IsingSampler(round(propBoot*Np), noDiag(sampleResult$graph), sampleResult$intercepts)
              
            } else if (model == "GGM") {
              g <- -sampleResult$graph
              diag(g) <- 1
              bootData <- mvtnorm::rmvnorm(round(propBoot*Np), sigma = corpcor::pseudoinverse(g))
              
            } else if (model == "graphicalVAR"){
              
              stop("model = 'graphicalVAR' not yet supported")
              
              
            } else stop(paste0("Model '",model,"' not supported."))
            
          } else {
            nPerson <- Np
            
            if (datatype == "normal"){
              bootData <- data[sample(seq_len(Np), round(propBoot*Np), replace=replacement), ]                      
            } else {
              bootData <- data
              bootSample <- sample(seq_len(Np), round(propBoot*Np), replace=replacement)
              bootData$data_c <- bootData$data_c[bootSample,]
              bootData$data_l <- bootData$data_l[bootSample,]
            }
            
          }
          
        } else if (type == "node") {
          
          # Nodewise
          nPerson <- Np
          nNode <- sample(subNodes,1)
          inSample <- sort(sample(seq_len(N),nNode))
          if (datatype == "normal"){
            bootData <- data[,inSample, drop=FALSE]            
          } else if (datatype == "graphicalVAR") {
            bootData <- data
            bootData$data_c <- bootData$data_c[,data$vars[inSample], drop = FALSE]
            bootData$data_l <- bootData$data_l[,c("1",grep(data$vars[inSample],names(data$data_l),value=TRUE)), drop = FALSE]
          }
          
        } else {
          # Personwise:
          nPerson <- sample(subCases,1)
          inSample <- 1:N
          persSample <- sort(sample(seq_len(Np),nPerson))
          if (datatype == "normal"){
            bootData <- data[persSample,, drop=FALSE]            
          } else if (datatype == "graphicalVAR") {
            bootData <- data
            bootData$data_c <- bootData$data_c[persSample,, drop = FALSE]
            bootData$data_l <- bootData$data_l[persSample,, drop = FALSE]
          }
        }
        
        # Some checks to remove progress bars:
        # if (!missing(prepFun)){
        #   # EBICglasso:
        #   if (!missing(prepArgs) & is.list(prepArgs) & identical(prepFun,qgraph::cor_auto)){
        #     prepArgs$verbose <- FALSE
        #   }
        # }
        # 
        res <- suppressWarnings(try({
          # estimateNetwork(bootData, 
          #                 default = default,
          #                 fun = fun,
          #                 prepFun = prepFun, # Fun to produce the correlation or covariance matrix
          #                 prepArgs = prepArgs, # list with arguments for the correlation function
          #                 estFun = estFun, # function that results in a network
          #                 estArgs = estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
          #                 graphFun = graphFun, # set to identity if missing
          #                 graphArgs = graphArgs, # Set to null if missing
          #                 intFun = intFun, # Set to null if missing
          #                 intArgs = intArgs, # Set to null if missing
          #                 labels = labels[inSample],
          #                 verbose = FALSE,
          #                 construct = construct)
          # 
          
          estimateNetwork(bootData, 
                          default = default,
                          fun = inputCheck$estimator,
                          .dots = inputCheck$arguments,
                          labels = labels[inSample],
                          verbose = FALSE,
                          weighted = weighted,
                          signed = signed,
                          .input = inputCheck,
                          memorysaver = memorysaver)
          
        }))
        if (is(res, "try-error")){
          
          if (tryCount == tryLimit) {
            stop("Maximum number of errors in bootstraps reached")
          }
          
          # warning("Error in bootstrap; retrying")
          tryCount <- tryCount + 1
        } else {
          break
        }
        
      }
      
      bootResults[[b]] <- res
      
      if (verbose){
        setTxtProgressBar(pb, b)
      }
    }
    if (verbose){
      close(pb)
    }
  } else {
    if (verbose){
      message("Bootstrapping...")
    }
    nClust <- nCores - 1
    cl <- makePSOCKcluster(nClust)
    
    # IF graph or data is missing, dummy graph:
    if (missing(graph)){
      graph <- matrix(0,N,N)
    }
    if (missing(data)){
      data <- matrix(0,Np,N)
    }
    if (missing(intercepts)){
      intercepts <- rep(0,N)
    }
    if (missing(sampleSize)){
      sampleSize <- Np
    }
    
    # Needed arguments to be excluded:
    excl <- c("prepFun", "prepArgs", "estFun", "estArgs", "graphFun", 
              "graphArgs", "intFun", "intArgs", "fun")
    
    clusterExport(cl, ls()[!ls()%in%c(excl,"cl")], envir = environment())
    # clusterExport(cl, export, envir = environment())
    
    
    # Run loop:
    bootResults <- pblapply(seq_len(nBoots), function(b){
      # Set library:
      .libPaths(library)
      
      tryLimit <- 10
      tryCount <- 0
      repeat{
       
        if (! type %in% c("node","person")){
          nNode <- N
          inSample <- seq_len(N)
          
          if (type == "jackknife"){
            if (datatype == "normal"){
              bootData <- data[-b,,drop=FALSE]                  
            } else {
              bootData <- data
              bootData$data_c <- bootData$data_c[-b,,drop=FALSE]
              bootData$data_l <- bootData$data_l[-b,,drop=FALSE]    
            }
            
            nPerson <- Np - 1
          } else if (type == "parametric"){
            nPerson <- Np
            if (model == "Ising"){
              bootData <- IsingSampler(round(propBoot*Np), noDiag(sampleResult$graph), sampleResult$intercepts)
              
            } else if (model == "GGM") {
              g <- -sampleResult$graph
              diag(g) <- 1
              bootData <- mvtnorm::rmvnorm(round(propBoot*Np), sigma = corpcor::pseudoinverse(g))
              
            } else if (model == "graphicalVAR"){
              
              stop("model = 'graphicalVAR' not yet supported")
              
              
            } else stop(paste0("Model '",model,"' not supported."))
            
          } else {
            nPerson <- Np
            
            if (datatype == "normal"){
              bootData <- data[sample(seq_len(Np), round(propBoot*Np), replace=replacement), ]                      
            } else {
              bootData <- data
              bootSample <- sample(seq_len(Np), round(propBoot*Np), replace=replacement)
              bootData$data_c <- bootData$data_c[bootSample,]
              bootData$data_l <- bootData$data_l[bootSample,]
            }
            
          }
          
        } else if (type == "node") {
          
          # Nodewise
          nPerson <- Np
          nNode <- sample(subNodes,1)
          inSample <- sort(sample(seq_len(N),nNode))
          if (datatype == "normal"){
            bootData <- data[,inSample, drop=FALSE]            
          } else if (datatype == "graphicalVAR") {
            bootData <- data
            bootData$data_c <- bootData$data_c[,data$vars[inSample], drop = FALSE]
            bootData$data_l <- bootData$data_l[,c("1",grep(data$vars[inSample],names(data$data_l),value=TRUE)), drop = FALSE]
          }
          
        } else {
          # Personwise:
          nPerson <- sample(subCases,1)
          inSample <- 1:N
          persSample <- sort(sample(seq_len(Np),nPerson))
          if (datatype == "normal"){
            bootData <- data[persSample,, drop=FALSE]            
          } else if (datatype == "graphicalVAR") {
            bootData <- data
            bootData$data_c <- bootData$data_c[persSample,, drop = FALSE]
            bootData$data_l <- bootData$data_l[persSample,, drop = FALSE]
          }
        }
        
        # Some checks to remove progress bars:
        # if (!missing(prepFun)){
        #   # EBICglasso:
        #   if (!missing(prepArgs) & is.list(prepArgs) & identical(prepFun,qgraph::cor_auto)){
        #     prepArgs$verbose <- FALSE
        #   }
        # }
        
        res <- suppressWarnings(try({
          # estimateNetwork(bootData,
          #                 default = default,
          #                 fun = fun,
          #                 prepFun = prepFun, # Fun to produce the correlation or covariance matrix
          #                 prepArgs = prepArgs, # list with arguments for the correlation function
          #                 estFun = estFun, # function that results in a network
          #                 estArgs = estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
          #                 graphFun = graphFun, # set to identity if missing
          #                 graphArgs = graphArgs, # Set to null if missing
          #                 intFun = intFun, # Set to null if missing
          #                 intArgs = intArgs, # Set to null if missing
          #                 labels = labels[inSample],
          #                 verbose = FALSE,
          #                 construct = construct)
          # 
          
          estimateNetwork(bootData, 
                          default = default,
                          fun = inputCheck$estimator,
                          .dots = inputCheck$arguments,
                          labels = labels[inSample],
                          verbose = FALSE,
                          weighted = weighted,
                          signed = signed,
                          .input = inputCheck,
                          memorysaver = memorysaver)
          
        }))
        if (is(res, "try-error")){
          
          if (tryCount == tryLimit) {
            stop("Maximum number of errors in bootstraps reached")
          }
          
          # warning("Error in bootstrap; retrying")
          tryCount <- tryCount + 1
        } else {
          break
        }
        
      }
      
      return(res)
    }, cl = cl)
  }
  
  
  #   if (edgeResample){
  #     bootGraphs <- do.call(abind::abind,c(lapply(bootResults,'[[','graph'),along=3))
  #     
  #     sampleDistribution <- sort(sampleGraph[upper.tri(sampleGraph,diag=FALSE)])
  #     for (b in seq_along(bootResults)){
  #       bootEdges <- bootResults[[b]]$graph[upper.tri(bootResults[[b]]$graph,diag=FALSE)]
  #       bootRank <- order(order(bootEdges))
  #       bootResults[[b]]$graph[upper.tri(bootResults[[b]]$graph,diag=FALSE)] <- sampleDistribution[bootRank]
  #       bootResults[[b]]$graph[lower.tri(bootResults[[b]]$graph,diag=FALSE)] <- t(bootResults[[b]]$graph)[lower.tri(bootResults[[b]]$graph,diag=FALSE)] 
  #     }
  #     
  #   }
  
  ### Compute the full parameter table!!
  if (verbose){
    message("Computing statistics...")
  }
  statTableOrig <- statTable(sampleResult,  name = "sample", alpha = alpha, computeCentrality = computeCentrality,statistics=statistics, directed=directed, includeDiagonal=includeDiagonal, communities=communities, useCommunities=useCommunities)
  
  if (nCores == 1){
    if (verbose){
      pb <- txtProgressBar(0,nBoots,style = 3)
    }
    statTableBoots <- vector("list", nBoots)
    for (b in seq_len(nBoots)){
      statTableBoots[[b]] <- statTable(bootResults[[b]], name = paste("boot",b), alpha = alpha, computeCentrality = computeCentrality, statistics=statistics, directed=directed,  communities=communities, useCommunities=useCommunities,includeDiagonal=includeDiagonal)
      if (verbose){
        setTxtProgressBar(pb, b)
      }
    }
    if (verbose){
      close(pb)
    }
  }  else {
    statTableBoots <- pblapply(seq_len(nBoots),function(b){
      # Set library:
      .libPaths(library)
      
      statTable(bootResults[[b]], name = paste("boot",b), alpha = alpha, computeCentrality = computeCentrality, statistics=statistics, directed=directed, communities=communities, useCommunities=useCommunities,includeDiagonal=includeDiagonal)
    }, cl = cl)
    # Stop the cluster:
    stopCluster(cl)
  }
  
  
  # Ordereing by node name to make nice paths:
  Result <- list(
    sampleTable = ungroup(statTableOrig),
    bootTable =  ungroup(dplyr::bind_rows(statTableBoots)),
    sample = sampleResult,
    boots = bootResults,
    type = type,
    sampleSize = Np)
  
  class(Result) <- "bootnet"
  
  return(Result)
  #   } else {
  #     
  #     ### Nodewise bootstrapping!
  #     
  #     # Bootstrap results:
  #     bootResults <- vector("list", nBoots)
  #     
  #     # Original centrality:
  #     origCentrality <- centrality(sampleResult$graph)
  #     
  #     # Setup the bootstrap table
  #     N <- ncol(data)
  #     simResults <- data.frame(id = seq_len(nBoots), nNodes = sample(nNodes,nBoots,TRUE))
  #     simResults[c("corStrength","corBetweenness","corCloseness","corSPL")] <- NA
  #     Strength <- Closeness <- Betweenness <- matrix(NA, nrow(simResults), N)
  #     colnames(Strength) <- colnames(Closeness) <- colnames(Betweenness) <- labels
  #     
  #     
  #     if (verbose){
  #       message("Bootstrapping...")
  #       pb <- txtProgressBar(0,nBoots,style = 3)
  #     }
  #     
  #     
  #     for (b in seq_len(nBoots)){
  #       nNodes <- simResults$nNodes[b]
  #       inSample <- sort(sample(seq_len(N),nNodes))
  #       bootData <- data[,inSample, drop=FALSE]
  #       res <- estimateNetwork(bootData, prepFun, prepArgs, estFun, estArgs)
  #       bootResults[[b]] <- list(
  #         graph = do.call(graphFun,c(list(res), graphArgs)),
  #         intercepts = do.call(intFun,c(list(res), intArgs)),
  #         results = res,
  #         labels = labels
  #       )
  #       
  #       class(bootResults[[b]]) <- c("bootnetResult", "list")
  #       
  #       if (verbose){
  #         setTxtProgressBar(pb, b) 
  #       }
  #     }
  #     
  #     if (verbose){
  #       close(pb)
  #     }
  #     
  #     browser()
  #     
  #   
  #       simCentrality <- centrality(bootResults[[b]]$graph)
  #       simResults$corStrength[b] <- cor(origCentrality$OutDegree[inSample], simCentrality$OutDegree)
  #       simResults$corBetweenness[b] <- cor(origCentrality$Betweenness[inSample], simCentrality$Betweenness)
  #       simResults$corCloseness[b] <- cor(origCentrality$Closeness[inSample], simCentrality$Closeness)
  #       simResults$corSPL[b] <- cor(origCentrality$ShortestPathLengths[inSample,inSample][upper.tri(origCentrality$ShortestPathLengths[inSample,inSample], diag=FALSE)], simCentrality$ShortestPathLengths[upper.tri(simCentrality$ShortestPathLengths, diag=FALSE)])
  #       
  #       Strength[b,inSample] <- simCentrality$OutDegree
  #       Closeness[b,inSample] <- simCentrality$Closeness
  #       Betweenness[b, inSample] <- simCentrality$Betweenness
  #       
  # 
  #     
  #     
  #     browser()
  #     
  #   }
}