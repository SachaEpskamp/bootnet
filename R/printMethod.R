getRefs <- function(x){
  citation <- switch(
    x,
    "none" = "",
    "EBICglasso" = c("Friedman, J. H., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9 (3), 432-441.",
                     "Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models. , 23 , 2020-2028.",
                     "Friedman, J. H., Hastie, T., & Tibshirani, R. (2014). glasso: Graphical lasso estimation of gaussian graphical models. Retrieved from https://CRAN.R-project.org/package=glasso",
                     "Epskamp, S., Cramer, A., Waldorp, L., Schmittmann, V. D., & Borsboom, D. (2012). qgraph: Network visualizations of relationships in psychometric data. Journal of Statistical Software, 48 (1), 1-18."
    ),
    "glasso" = c("Friedman, J. H., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9 (3), 432-441.",
                 "Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models. , 23 , 2020-2028.",
                 "Friedman, J. H., Hastie, T., & Tibshirani, R. (2014). glasso: Graphical lasso estimation of gaussian graphical models. Retrieved from https://CRAN.R-project.org/package=glasso",
                 "Epskamp, S., Cramer, A., Waldorp, L., Schmittmann, V. D., & Borsboom, D. (2012). qgraph: Network visualizations of relationships in psychometric data. Journal of Statistical Software, 48 (1), 1-18."
    ),
    "IsingFit" = "van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo, L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for constructing networks from binary data. Scientific reports, 4 (5918), 1-10.",
    "IsingLL" = c("Epskamp, S., Maris, G., Waldorp, L., & Borsboom, D. (in press). Network psychometrics. In P. Irwing, D. Hughes, & T. Booth (Eds.), Handbook of psychometrics. New York, NY, USA: Wiley.",
                  "Epskamp, S. (2014). IsingSampler: Sampling methods and distribution functions for the Ising model. Retrieved from github.com/SachaEpskamp/IsingSampler"),
    "huge" = "Zhao, T., Li, X., Liu, H., Roeder, K., Lafferty, J., & Wasserman, L. (2015). huge: High-dimensional undirected graph estimation. Retrieved from https://CRAN.R-project.org/package=huge",
    "adalasso" = "Kraeamer, N., Schaeafer, J., & Boulesteix, A.-L. (2009). Regularized estimation of large-scale gene association networks using graphical gaussian models. BMC Bioinformatics, 10 (1), 1-24."
  )
  
  citation <- c(citation,
                "Epskamp, S., Borsboom, D., & Fried, E. I. (2016). Estimating psychological networks and their accuracy: a tutorial paper. arXiv preprint, arXiv:1604.08462.")
  
  citation
}

print.bootnet <- function(x, ...){

  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  cat("=== bootnet Results ===")
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
      paste0("\nUse plot(",name,") to plot summarized statistics (see ?plot.bootnet for details)"),
      "\n\nRelevant references:\n\n",paste0("\t",getRefs(x$sample$input$default),collapse="\n")
      )
}

print.bootnetResult <- function(x, ...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  cat("=== Estimated network ===")
  cat("\nNumber of nodes:",nrow(x[['graph']]),
      "\nNumber of non-zero edges:",sum(x[['graph']][upper.tri(x[['graph']],diag=FALSE)]==0) ,
      "\nSparsity:",mean(x[['graph']][upper.tri(x[['graph']],diag=FALSE)]==0) ,
      paste0("\nNetwork stored in ",name,"$graph"),
      "\n",
      paste0("\nDefault set used: ",x$input$default),     
      "\n",
      paste0("\nUse plot(",name,", layout = 'spring') to plot estimated network"),
      paste0("\nUse bootnet(",name,") to bootstrap edge weights and centrality indices"),
      "\n\nRelevant references:\n\n",paste0("\t",getRefs(x$input$default),collapse="\n")
  )
}