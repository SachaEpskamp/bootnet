# Function to run a standard network simulation study:
# Two methods: dataGenerator function, or using graph, intercepts and model for GGM/Ising.

replicationSimulator <- function(
  input = genGGM(Nvar = 10),  # A matrix, or a list with graph and intercepts elements.
  nCases = c(50,100,250,500,1000,2500), # Number of cases
  nReps = 100, # Number of repititions per condition
  nCores = 1, # Number of computer cores used
  default,
  dataGenerator, 
  ..., # estimateNetwork arguments (if none specified, will default to default = "EBICglasso)
  moreArgs = list() # List of extra args not intended to be varied as conditions
){
  # Dots list:
  .dots <- list(...)
  # default <- match.arg(default)
  
  # Check default and dataGenerator:
  if (missing(default) & missing(dataGenerator)){
    message("'default' and 'dataGenerator' are missing. Setting default = 'EBICglasso'")
    default <- "EBICglasso"
  }
  
  if (missing(default) & length(.dots) == 0){
    message("No estimator specified. Setting default = 'EBICglasso'")
    default <- "EBICglasso"
  }
  
  # Data generator:
  if (!missing(default) & missing(dataGenerator)){
    
    if (default == "EBICglasso" || default == "glasso" || default == "pcor" || default == "adalasso" || default == "huge" || default == "ggmModSelect"){
      message("Setting 'dataGenerator = ggmGenerator(ordinal = FALSE)'")
      dataGenerator <- ggmGenerator(ordinal = FALSE)
    } else if (default == "IsingFit" || default == "IsingSampler"){
      message("Setting 'dataGenerator = IsingGenerator()'")
      dataGenerator <- IsingGenerator()
    } else {
      stop(paste0("Default set '",default, "' not yet supported. Please manually specify 'dataGenerator'"))
    }
    
    
  }
  
  # parSim arguments:
  Args <- c(
    list(
      # Conditions:
      nCases = nCases,
      default = default,
      
      # Setup:,
      write=FALSE,
      nCores = nCores,
      reps = nReps,
      debug=FALSE,
      export = c("input","dataGenerator",".dots","moreArgs"),
      
      expression = expression({
        cor0 <- function(x,y){
          if (all(is.na(x)) || all(is.na(y))){
            return(NA)
          }
          
          if (sd(x,na.rm=TRUE)==0 | sd(y,na.rm=TRUE) == 0){
            return(0)
          }
          
          return(cor(x,y,use="pairwise.complete.obs"))
        }
        
        # Generate the input:
        if (is.function(input)){
          inputResults <- input()
        } else {
          inputResults <- input
        }
        
        # True network:
        if (is.list(inputResults)){
          trueNet <- inputResults$graph
        } else {
          trueNet <- inputResults
        }
        
        # Generate the datasets:
        Data1 <- dataGenerator(nCases, inputResults)
        Data2 <- dataGenerator(nCases, inputResults)
        
        # Compute network 1:
        args1 <- list()
        args1$data <- Data1
        args1$verbose <- FALSE
        args1$default <- default
        args1 <- c(args1,moreArgs)
        
        for (i in seq_along(.dots)){
          args1[[names(.dots)[i]]] <- get(names(.dots)[i])
        }
        netResults1 <- do.call(bootnet::estimateNetwork,args1)
        estNet1 <- qgraph::getWmat(netResults1)
        
        
        # Compute network 2:
        args2 <- list()
        args2$data <- Data2
        args2$verbose <- FALSE
        args2$default <- default
        args2 <- c(args2,moreArgs)
        
        for (i in seq_along(.dots)){
          args2[[names(.dots)[i]]] <- get(names(.dots)[i])
        }
        netResults2 <- do.call(bootnet::estimateNetwork,args2)
        estNet2 <- qgraph::getWmat(netResults2)
        
        
        # Compute measures:
        ### STORE RESULTS ###
        SimulationResults <- list()

        # Estimated edges:
        est1 <- estNet1[upper.tri(estNet1)]
        est2 <- estNet2[upper.tri(estNet2)]
        
        # Equal?
        SimulationResults$identical <- all((est1 == 0) == (est2 == 0))
        
        # Correlation between edge wegights:
        SimulationResults$correlation <- cor(est1, est2)
        
        # Correlation nonzero edge wegights:
        SimulationResults$correlationNonZero <- cor(est1[est1!=0 & est2 != 0], est2[est1!=0 & est2 != 0])
        
        # Jaccard index:
        SimulationResults$jaccard <-  sum(est1!=0 & est2!=0) / sum(est1!=0 | est2!=0) 
        
        # Percentage of edges in network 1 replicated in network 2:
        SimulationResults$replicatedEdges <- sum(est1!=0 & est2!=0)/sum(est1!=0)
        
        # Percentage of zeroes in network 1 replicated in network 2:
        SimulationResults$replicatedZeroes <- sum(est1==0 & est2==0)/sum(est1==0)
        
        # Centrality:
        centEst1 <- qgraph::centrality(estNet1)
        centEst2 <- qgraph::centrality(estNet2)
        
        SimulationResults$strength <- cor0(centEst1$OutDegree,centEst2$OutDegree)
        SimulationResults$closeness <- cor0(centEst1$Closeness,centEst2$Closeness)
        SimulationResults$betweenness <- cor0(centEst1$Betweenness,centEst2$Betweenness)
   
        SimulationResults
      })),
      .dots)

    Results <- do.call(parSim, Args)
    
    class(Results) <- c("replicationSimulator","data.frame")
    return(Results)
}

