cor0 <- function(x,y,...){
  if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2 || sd(x,na.rm=TRUE)==0 | sd(y,na.rm=TRUE) == 0){
    return(0)
  } else {
    return(cor(x,y,...))
  }
}

# Smetric
corStability <- function(x,cor=0.7, statistics = c("strength","closeness","betweenness"),
                         verbose = TRUE){

  if (!x$type %in% c("node","person")){
    stop("S-metric only available for person or node drop bootstrap")
  }
  
  if (x$type == "node"){
    x$bootTable$prop <- 1 - (x$bootTable$nNode / x$sample$nNodes)
  } else {
    x$bootTable$prop <- 1 - (x$bootTable$nPerson / x$sample$nPerson)
  }
  
  sample <- x$sampleTable %>% 
    filter(type %in% statistics) %>%
    select(node1,node2,type,original = value)
    
  
  
  max0 <- function(x){
    if (length(x)==0)return(0) else return(max(x))
  }

  S <- x$bootTable %>%
    filter(type %in% statistics) %>%
    left_join(sample,by=c("node1","node2","type")) %>% 
    group_by(name,type,prop) %>% 
    summarize(stability = cor0(value,original)) %>% 
    group_by(prop,type) %>% 
    summarize(P = mean(stability > cor)) %>% 
    group_by(type) %>%
    summarize(Smetric = max0(prop[P > 0.95]))
  
  # all unique sampling levels:
  samplingLevels <- sort(unique(x$bootTable$prop))
  
  Smetric <- S$Smetric
  names(Smetric) <- S$type
  
  # Print information per sampling level:
  if (verbose){
    # Get counts:
    counts <- x$bootTable %>% 
      group_by(name) %>% summarize(nPerson = unique(nPerson)) %>%
      group_by(nPerson) %>%
      tally %>% arrange(nPerson) %>% 
      as.data.frame
    
    counts[['Drop%']] <- round(100 * (1 - counts$nPerson / x$sample$nPerson),1)
    rownames(counts) <- NULL
    counts <- counts[,c("nPerson","Drop%","n")]
    
    cat("=== Correlation Stability Analysis ===",
        "\n\nSampling levels tested:\n")
    print(counts)
    cat(paste0("\nMaximum drop proportions to retain correlation of ",cor," in at least 95% of the samples:\n\n"))
    
    samplingLevels <- c(0,samplingLevels,1)
    
    
    for (i in seq_along(Smetric)){
      if (any(samplingLevels < Smetric[i])){
        lower <- max(which(samplingLevels < Smetric[i]))
      } else {
        lower <- 1
      }
      
      if (any(samplingLevels > Smetric[i])){
        upper <- min(which(samplingLevels > Smetric[i]))
      } else {
        upper <- length(samplingLevels)
      }
        
      
      cat(paste0(names(Smetric)[i],": ",round(Smetric[i],3),
          "\n  - For more accuracy, run bootnet(..., caseMin = ",round(samplingLevels,3)[lower],", caseMax = ",round(samplingLevels,3)[upper],")"),
          "\n\n"
      )
      
    }
    cat("Accuracy can also be increased by increasing both 'nBoots' and 'caseN'.")
    
  }
  
  
  invisible(Smetric)
}