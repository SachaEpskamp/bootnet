# Null function:
null <- function(...) NULL

### IsingFit ###
# prep fun (make binary if needed):
binarize <- function(x, split = "median", na.rm=TRUE, removeNArows = TRUE){
  x <- as.data.frame(x)
  
  if (all(unlist(x) %in% c(0,1))){
    return(x)
  } else {
    
    if (is.function(split)){
      splitName <- deparse(substitute(split))
    } else {
      splitName <- split
    }
    warning(paste("Splitting data by",splitName))
     
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