print.netSimulator <- function(x,...) summary(x,...)

summary.netSimulator <- function(object, ...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  
  Exclude <- c(
    "rep","id","correctModel","sensitivity","specificity","correlation","error","errorMessage"
  )
  # check number of levels:
  Conditions <- names(object)[!names(object)%in%Exclude]
  
  . <- NULL
  
  # Summarize per case:
  df <- object %>% dplyr::select_("sensitivity","specificity","correlation",.dots = Conditions) %>% 
    dplyr::group_by_(.dots = Conditions) %>% dplyr::summarize_each(funs(mean(.,na.rm=TRUE))) %>% 
    dplyr::arrange_(~nCases) %>% as.data.frame
  
  dfSD <- object %>% dplyr::select_("sensitivity","specificity","correlation",.dots = Conditions) %>% 
    dplyr::group_by_(.dots = Conditions) %>% dplyr::summarize_each(funs(sd(.,na.rm=TRUE))) %>% 
    dplyr::arrange_(~nCases) %>% as.data.frame
  
  row.names(df) <- NULL
  
  cat("=== netSimulator Results ===\n\n")
  cat("Mean values per varied levels:\n\n")
  print(df)
  cat("Standard deviation per varied levels:\n\n")
  print(dfSD)
  
  
  cat(paste0("\n\nUse plot(",name,") to plot results (nCases only), or as.data.frame(",name,") to see all results."))
  invisible(df)
}

plot.netSimulator <- function(x, print = TRUE, ylim = c(0,1), ...){
  # Gather:
  Gathered <- x %>% dplyr::select_("nCases","sensitivity","specificity","correlation") %>%
    gather_("index","value",c("sensitivity","specificity","correlation"))
  
  # Plot:
  g <- ggplot(Gathered, aes_string(x="factor(nCases)",y="value")) + 
    facet_grid(~ index) + geom_boxplot() + theme_bw() + ylim(ylim[1],ylim[2]) + 
    ylab("") + xlab("Number of cases")
  
  if (print){
    print(g)
    invisible(g)
  } else {
    return(g)
  }
}