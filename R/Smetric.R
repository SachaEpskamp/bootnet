cor0 <- function(x,y,...){
  if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2 || sd(x,na.rm=TRUE)==0 | sd(y,na.rm=TRUE) == 0){
    return(0)
  } else {
    return(cor(x,y,...))
  }
}

# Smetric
corStability <- function(x,cor=0.7, statistics = c("strength","closeness","betweenness")){
  
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
  
  Smetric <- S$Smetric
  names(Smetric) <- S$type
  return(Smetric)
}