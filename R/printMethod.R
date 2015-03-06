print.bootnet <- function(x, ...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  cat("=== bootnet results ===")
  cat("\nNumber of nodes:",nrow(x[['sample']][['graph']]),
      "\nNumber of non-zero edges in sample:",sum(x[['sample']][['graph']][upper.tri(x[['sample']][['graph']],diag=FALSE)]==0) ,
      "\nSparsity in sample:",mean(x[['sample']][['graph']][upper.tri(x[['sample']][['graph']],diag=FALSE)]) ,
      "\nNumber of intercepts:",NROW(x[['sample']][['intercepts']]),
      "\nNumber of bootstrapped networks:",length(x[['boots']]),
      paste0("\nResults of original sample stored in ",name,"$sample"),
      paste0("\nTable of all statistics from original sample stored in ",name,"$sampleTable"),
      paste0("\nResults of bootstraps stored in ",name,"$boots"),
      paste0("\nTable of all statistics from bootstraps stored in ",name,"$bootTable"),
      "\n",
      paste0("\nUse plot(",name,"$sample, layout = 'spring') to plot estimated network of original sample"),
      paste0("\nUse summary(",name,") to inspect summarized statistics (see ?summary.bootnet for details)"),
      paste0("\nUse plot(",name,") to plot summarized statistics (see ?plot.bootnet for details)")
      )
}