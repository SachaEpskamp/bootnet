# Copy-paste from psychonetrics:
quantiletransform <- function(x){
  xNoNA <- x[!is.na(x)]
  ord <- order(xNoNA)
  sorted <- sort(xNoNA)
  nBelow <- rank(sorted, ties.method = "min")
  p <- nBelow / (max(nBelow)+1)
  q <- qnorm(p)
  xTrans <- x
  xTrans[!is.na(xTrans)][ord] <- q
  return(xTrans)
}

quantile_transformation <- function(x){
  for (i in 1:ncol(x)){
    x[,i] <- quantiletransform(x[,i])
  }
  x
}

rank_transformation <- function(x,ties.method = c("average", "first", "last", "random", "max", "min")){
  ties.method <- match.arg(ties.method)
  for (i in 1:ncol(x)){
    x[!is.na(x[,i]),i] <- rank(x[!is.na(x[,i]),i], ties.method = ties.method)
  }
  x
}