# Function to threshold network based on bootstrap samples:
bootThreshold <- function(bootobject, alpha = 0.05,verbose=TRUE){
  # Check if object is bootnet object:
  if (class(bootobject) != "bootnet"){
    stop("'bootobject' must be an object of class 'bootnet'")
  }
  
  # Check type:
  if (bootobject$type != "nonparametric" & bootobject$type != "parametric"){
    stop("Bootstrap type must be 'nonparametric' or 'parametric'")
  }
  
  # Check alpha:
  if (verbose){
    exp <- expAlpha(alpha,length(bootobject$boots))
    message(paste0("Expected significance level given number of bootstrap samples is approximately: ",format(signif(exp,2),scientific = FALSE)))
  }
  
  # Extract the network object:
  Network <- bootobject$sample
  
  # Summary table of edge weights:
  bootSummary <- bootobject$bootTable %>% 
    dplyr::filter_(~type == "edge") %>%
    dplyr::group_by_(~node1,~node2) %>%
    dplyr::summarize(
      lower = quantile(value, alpha/2),
      upper = quantile(value, 1 - alpha/2)
    ) %>% 
    dplyr::mutate(sig = upper < 0 | lower > 0) %>%
    filter_(~!sig)
  
  # Threshold network:
  if (nrow(bootSummary) > 0){
    for (i in 1:nrow(bootSummary)){
      Network$graph[Network$labels == bootSummary$node1[i],Network$labels == bootSummary$node2[i]] <- 0
      if (!Network$directed){
        Network$graph[Network$labels == bootSummary$node2[i],Network$labels == bootSummary$node1[i]] <- 0
      }
    }
  } else {
    message("All edges indicated to be nonzero")
  }
  
  # Add indicator network is thresholded:
  Network$thresholded <- TRUE
  
  # Return network:
  return(Network)
}