print.replicationSimulator <- function(x, digits = 2, ...) summary(x, digits = digits, ...)

summary.replicationSimulator <- function(object, digits = 2, ...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  
  # Check for errors:
  if (all(object$error)) stop(paste0("All simulations resulted in errors:\n",paste(unique(object$errorMessage, collapse = "\n"))))
  
  Exclude <- c("rep", "id", "identical", "correlation", 
               "correlationNonZero", "jaccard", "replicatedEdges", "replicatedZeroes", 
               "strength", "closeness", "betweenness", "error", "errorMessage")
 
  # check number of levels:
  Conditions <- names(object)[!names(object)%in%Exclude]
  
  . <- NULL
  
  fun <- function(x,digits=2){
    paste0(round(mean(x,na.rm=TRUE),digits), " (",round(sd(x,na.rm=TRUE),digits),")")
  }
  
  # Summarize per case:
  df <- object %>% dplyr::select_("correlation", 
                                  "correlationNonZero", "jaccard", "replicatedEdges", "replicatedZeroes", 
                                  "strength", "closeness", "betweenness",.dots = Conditions) %>% 
    dplyr::group_by_(.dots = Conditions) %>% dplyr::summarize_each(funs(fun(.,digits=digits))) %>% 
    dplyr::arrange_(~nCases) %>% as.data.frame
  # 
  # dfSD <- object %>% dplyr::select_("sensitivity","specificity","correlation",.dots = Conditions) %>% 
  #   dplyr::group_by_(.dots = Conditions) %>% dplyr::summarize_each(funs(sd(.,na.rm=TRUE))) %>% 
  #   dplyr::arrange_(~nCases) %>% as.data.frame
  
  row.names(df) <- NULL
  
  cat("=== replicationSimulator Results ===\n\n")
  cat("Mean (SD) values per varied levels:\n\n")
  print(df)
  # cat("Standard deviation per varied levels:\n\n")
  # print(dfSD)
  
  
  cat(paste0("\n\nUse plot(",name,") to plot results (nCases only), or as.data.frame(",name,") to see all results."))
  invisible(df)
}

# Plot method

plot.replicationSimulator <- function(x, yvar = c("correlation","jaccard","replicatedEdges","replicatedZeroes"), ...){
  plot.netSimulator(x, yvar = yvar, ...)
}