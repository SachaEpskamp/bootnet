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
  default = c("none", "EBICglasso", "pcor","IsingFit","IsingLL"), # Default method to use. EBICglasso, IsingFit, concentration, some more....
  type = c("nonparametric","parametric","node","person","jackknife"), # Bootstrap method to use
  model = c("detect","GGM","Ising"), # Models to use for bootstrap method = parametric. Detect will use the default set and estimation function.
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
  alpha = 1, # centrality alpha
  subNodes = 2:(ncol(data)-1), # if type = "node", defaults to 2:(p-1)
  subPersons = round(seq(0.25,0.95,length=10) * nrow(data)),
  computeCentrality = TRUE,
  propBoot = 1, # M out of N
  # subsampleSize,
  replacement = TRUE,
  edgeResample = FALSE # If true, only resample edges from original estimate
  # scaleAdjust = FALSE
){
  if (default[[1]]=="glasso") default <- "EBICglasso"
  default <- match.arg(default)
  type <- match.arg(type)
  model <- match.arg(model)
  N <- ncol(data)
  Np <- nrow(data)
  
  if (type == "jackknife"){
    message("Jacknife overwrites nBoot to sample size")
    nBoots <- Np
  }
  
  
  if (type == "node" & N < 3){
    stop("Node-wise bootstrapping requires at least three nodes.")
  }
  
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
                        pcor = qgraph::cor_auto
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
                         EBICglasso = ifElse(identical(prepFun,qgraph::cor_auto),list(verbose=FALSE),
                                             ifelse(identical(prepFun,cor),list(use = "pairwise.complete.obs"),list())),
                         IsingFit = list(),
                         pcor =  ifElse(identical(prepFun,qgraph::cor_auto),list(verbose=FALSE),
                                        ifelse(identical(prepFun,cor),list(use = "pairwise.complete.obs"),list())),
                         IsingLL = list()
      )
      
      
    }
    
    # estFun:
    if (missing(estFun)){
      estFun <- switch(default,
                       EBICglasso = qgraph::EBICglasso,
                       pcor = corpcor::cor2pcor,
                       IsingFit = IsingFit::IsingFit,
                       IsingLL = IsingSampler::EstimateIsing
      )
    }
    
    # estArgs:
    if (missing(estArgs)){
      estArgs <- switch(default,
                        EBICglasso = list(n = Np, returnAllResults = TRUE),
                        IsingFit = list(plot = FALSE, progress = FALSE),
                        pcor = list(),
                        IsingLL = list(method = "ll")
      )
    }
    
    # graphFun:
    if (missing(graphFun)){
      graphFun <- switch(default,
                         EBICglasso = function(x)x[['optnet']],
                         IsingFit = function(x)x[['weiadj']],
                         pcor = function(x)as.matrix(Matrix::forceSymmetric(x)),
                         IsingLL = function(x)x[['graph']]
      )
    }
    
    # graphArgs:
    if (missing(graphArgs)){
      graphArgs <- switch(default,
                          EBICglasso = list(),
                          IsingFit = list(),
                          pcor = list(),
                          IsingLL = list()
      )
    }
    
    # intFun:
    if (missing(intFun)){
      intFun <- switch(default,
                       EBICglasso = null,
                       IsingFit = function(x)x[['thresholds']],
                       pcor = null,
                       IsingLL = function(x) x[['thresholds']]
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
  
  ## For parametric bootstrap, detect model
  if (type == "parametric" & model == "detect"){
    if (default != "none"){
      model <- ifelse(grepl("ising",default,ignore.case=TRUE),"Ising","GGM")
    } else {
      model <- ifelse(any(grepl("ising",deparse(estFun),ignore.case=TRUE)),"Ising","GGM")
    }
    message(paste0("model set to '",model,"'"))
  }
  
  # Estimate sample result:
  if (verbose){
    message("Estimating sample network...")
  }
  res <- estimateNetwork(data, prepFun, prepArgs, estFun, estArgs)
  sampleGraph <- do.call(graphFun,c(list(res), graphArgs))
  sampleResult <- list(
    graph = sampleGraph,
    intercepts = do.call(intFun,c(list(res), intArgs)),
    results = res,
    labels = labels,
    nNodes = ncol(data),
    nPerson = Np
  )
  class(sampleResult) <- c("bootnetResult", "list")
  
  if (!isSymmetric(sampleResult[['graph']])){
    stop("bootnet does not support directed graphs")
  }
  
  
  #   ### Observation-wise bootstrapping!
  #   if (type == "observation"){
  # Bootstrap results:
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
        nNode <- ncol(data)
        inSample <- seq_len(N)
        
        if (type == "jackknife"){
          bootData <- data[-b,,drop=FALSE]    
          nPerson <- Np - 1
        } else if (type == "parametric"){
          nPerson <- Np
          if (model == "Ising"){
            bootData <- IsingSampler(round(propBoot*Np), noDiag(sampleResult$graph), sampleResult$intercepts)
            
          } else if (model == "GGM") {
            g <- -sampleResult$graph
            diag(g) <- 1
            bootData <- mvtnorm::rmvnorm(round(propBoot*Np), sigma = corpcor::pseudoinverse(g))
            
          } else stop(paste0("Model '",model,"' not supported."))
          
        } else {
          nPerson <- Np
          bootData <- data[sample(seq_len(Np), round(propBoot*Np), replace=replacement), ]        
        }
        
      } else if (type == "node") {
        # Nodewise
        nPerson <- Np
        nNode <- sample(subNodes,1)
        inSample <- sort(sample(seq_len(N),nNode))
        bootData <- data[,inSample, drop=FALSE]
      } else {
        # Personwise:
        nNode <- ncol(data)
        nPerson <- sample(subPersons,1)
        inSample <- 1:N
        persSample <- sort(sample(seq_len(Np),nPerson))
        bootData <- data[persSample,,drop=FALSE]
         
      }
      
      res <- suppressWarnings(try(estimateNetwork(bootData, prepFun, prepArgs, estFun, estArgs)))
      if (is(res, "try-error")){
        if (tryCount == tryLimit) stop("Maximum number of errors in bootstraps reached")
        
        # warning("Error in bootstrap; retrying")
        tryCount <- tryCount + 1
      } else {
        break
      }
      
    }
    BootGraph <- do.call(graphFun,c(list(res), graphArgs))
    #     if (scaleAdjust == 1){
    #       if (!all(BootGraph[upper.tri(BootGraph,diag=FALSE)] == 0)){
    #         BootGraph[upper.tri(BootGraph,diag=FALSE)] <- BootGraph[upper.tri(BootGraph,diag=FALSE)] * sum(abs(sampleGraph[upper.tri(sampleGraph,diag=FALSE)])) /
    #           sum(abs(BootGraph[upper.tri(BootGraph,diag=FALSE)]))
    #         BootGraph[lower.tri(BootGraph,diag=FALSE)] <- t(BootGraph)[lower.tri(BootGraph,diag=FALSE)]        
    #       }
    #     }
    bootResults[[b]] <- list(
      graph = BootGraph,
      intercepts = do.call(intFun,c(list(res), intArgs)),
      results = res,
      labels = labels[inSample],
      nNodes = nNode,
      nPerson = nPerson
    )
    
    class(bootResults[[b]]) <- c("bootnetResult", "list")
    
    if (verbose){
      setTxtProgressBar(pb, b)
    }
  }
  if (verbose){
    close(pb)
  }
  
  
  #   if (isTRUE(scaleAdjust) || scaleAdjust == 2){
  #     if (verbose){
  #       message("Applying bias-adjustment correction")
  #     }
  #     bootGraphs <- do.call(abind::abind,c(lapply(bootResults,'[[','graph'),along=3))
  #     needAdjust <- apply(bootGraphs,1:2,function(x)min(x)<=0&max(x)>=0)
  #     needAdjust <- needAdjust[upper.tri(needAdjust,diag=FALSE)]
  #     if (any(needAdjust)){
  #       # adjust <- which(needAdjust & upper.tri(needAdjust,diag=FALSE),arr.ind=TRUE)
  #       sampleDistribution <- sort(sampleGraph[upper.tri(sampleGraph,diag=FALSE)])
  #       for (b in seq_along(bootResults)){
  #           bootEdges <- bootResults[[b]]$graph[upper.tri(bootResults[[b]]$graph,diag=FALSE)]
  #           bootRank <- order(order(bootEdges))
  #           bootResults[[b]]$graph[upper.tri(bootResults[[b]]$graph,diag=FALSE)] <- ifelse(needAdjust,sampleDistribution[bootRank],bootEdges)
  #           bootResults[[b]]$graph[lower.tri(bootResults[[b]]$graph,diag=FALSE)] <- t(bootResults[[b]]$graph)[lower.tri(bootResults[[b]]$graph,diag=FALSE)] 
  #       }
  #     }
  #     
  #   }
  if (edgeResample){
    bootGraphs <- do.call(abind::abind,c(lapply(bootResults,'[[','graph'),along=3))
    
    sampleDistribution <- sort(sampleGraph[upper.tri(sampleGraph,diag=FALSE)])
    for (b in seq_along(bootResults)){
      bootEdges <- bootResults[[b]]$graph[upper.tri(bootResults[[b]]$graph,diag=FALSE)]
      bootRank <- order(order(bootEdges))
      bootResults[[b]]$graph[upper.tri(bootResults[[b]]$graph,diag=FALSE)] <- sampleDistribution[bootRank]
      bootResults[[b]]$graph[lower.tri(bootResults[[b]]$graph,diag=FALSE)] <- t(bootResults[[b]]$graph)[lower.tri(bootResults[[b]]$graph,diag=FALSE)] 
    }
    
  }
  
  ### Compute the full parameter table!!
  if (verbose){
    message("Computing statistics...")
    pb <- txtProgressBar(0,nBoots+1,style = 3)
  }
  statTableOrig <- statTable(sampleResult,  name = "sample", alpha = alpha, computeCentrality = computeCentrality)
  if (verbose){
    setTxtProgressBar(pb, 1)
  }
  statTableBoots <- vector("list", nBoots)
  for (b in seq_len(nBoots)){
    statTableBoots[[b]] <- statTable(bootResults[[b]], name = paste("boot",b), alpha = alpha, computeCentrality = computeCentrality)    
    if (verbose){
      setTxtProgressBar(pb, b+1)
    }
  }
  if (verbose){
    close(pb)
  }
  
  # Ordereing by node name to make nice paths:
  Result <- list(
    sampleTable = ungroup(statTableOrig),
    bootTable =  ungroup(dplyr::rbind_all(statTableBoots)),
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