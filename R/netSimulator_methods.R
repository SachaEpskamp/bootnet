print.netSimulator <- function(x) summary(x)

summary.netSimulator <- function(object, ...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  
  # Summarize per case:
  df <- object %>% select(nCases,sensitivity,specificity,correlation) %>% 
    group_by(nCases) %>% summarize_each(funs(mean(.,na.rm=TRUE))) %>% 
    arrange(nCases) %>% as.data.frame
  
  row.names(df) <- NULL
  
  cat("=== netSimulator Results ===\n\n")
  cat("Mean values per sampling level:\n\n")
  print(df)
  cat(paste0("\n\nUse plot(",name,") to plot results."))
  invisible(df)
}

plot.netSimulator <- function(x, print = TRUE, ylim = c(0,1)){
  
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