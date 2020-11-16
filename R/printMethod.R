getRefs <- function(x){
  citation <- switch(
    x,
    "none" = "",
    "EBICglasso" = c("Friedman, J. H., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9 (3), 432-441.",
                     "Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models. ",
                     "Friedman, J. H., Hastie, T., & Tibshirani, R. (2014). glasso: Graphical lasso estimation of gaussian graphical models. Retrieved from https://CRAN.R-project.org/package=glasso",
                     "Epskamp, S., Cramer, A., Waldorp, L., Schmittmann, V. D., & Borsboom, D. (2012). qgraph: Network visualizations of relationships in psychometric data. Journal of Statistical Software, 48 (1), 1-18."
    ),
    "glasso" = c("Friedman, J. H., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9 (3), 432-441.",
                 "Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models. , 23 , 2020-2028.",
                 "Friedman, J. H., Hastie, T., & Tibshirani, R. (2014). glasso: Graphical lasso estimation of gaussian graphical models. Retrieved from https://CRAN.R-project.org/package=glasso",
                 "Epskamp, S., Cramer, A., Waldorp, L., Schmittmann, V. D., & Borsboom, D. (2012). qgraph: Network visualizations of relationships in psychometric data. Journal of Statistical Software, 48 (1), 1-18."
    ),
    "IsingFit" = "van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo, L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for constructing networks from binary data. Scientific reports, 4 (5918), 1-10.",
    "IsingSampler" = c("Epskamp, S., Maris, G., Waldorp, L., & Borsboom, D. (in press). Network psychometrics. In P. Irwing, D. Hughes, & T. Booth (Eds.), Handbook of psychometrics. New York, NY, USA: Wiley.",
                  "Epskamp, S. (2014). IsingSampler: Sampling methods and distribution functions for the Ising model. Retrieved from github.com/SachaEpskamp/IsingSampler"),
    "huge" = "Zhao, T., Li, X., Liu, H., Roeder, K., Lafferty, J., & Wasserman, L. (2015). huge: High-dimensional undirected graph estimation. Retrieved from https://CRAN.R-project.org/package=huge",
    "adalasso" = "Kraeamer, N., Schaeafer, J., & Boulesteix, A.-L. (2009). Regularized estimation of large-scale gene association networks using graphical gaussian models. BMC Bioinformatics, 10 (1), 1-24.",
    "mgm" = "Jonas M. B. Haslbeck, Lourens J. Waldorp (2016). mgm: Structure Estimation for Time-Varying Mixed Graphical Models in high-dimensional Data arXiv preprint:1510.06871v2 URL http://arxiv.org/abs/1510.06871v2.",
    "TMFG" = c("Christensen, A. P., Kenett, Y. N., Aste, T., Silvia, P. J., & Kwapil, T. R. (2018). Network Structure of the Wisconsin Schizotypy Scales-Short Forms: Examining Psychometric Network Filtering Approaches. Behavorial Research Methods. DOI: 10.3758/s13428-018-1032-9",
               "Christensen, A. P. (2018). NetworkToolbox: Methods and Measures for Brain, Cognitive, and Psychometric Network Analysis in R"),
    "LoGo" = c("Barfuss, W., Massara, G. P., Di Matteo, T., & Aste, T. (2016). Parsimonious modeling with information filtering networks. Physical Review E, 94(6), 062306.",
               "Christensen, A. P. (2018). NetworkToolbox: Methods and Measures for Brain, Cognitive, and Psychometric Network Analysis in R"),
    "ggmModSelect" = c("Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models.",
                       "http://psychosystems.org/qgraph_1.5"),
    "graphicalVAR" = c("Abegaz, F., & Wit, E. (2013). Sparse time series chain graphical models for reconstructing genetic networks. Biostatistics, 14(3), 586???599.",
                       "Epskamp, S., Waldorp, L. J., Mottus, R., & Borsboom, D. (2018). The Gaussian Graphical Model in Cross-sectional and Time-series Data. Multivariate Behavioral Research",
                       "Rothman, A. J., Levina, E., & Zhu, J. (2010). Sparse multivariate regression with covariance estimation. Journal of Computational and Graphical Statistics, 19(4), 947???962.",
                       "Wild, B., Eichler, M., Friederich, H.-C., Hartmann, M., Zipfel, S., & Herzog, W. (2010). A graphical vector autoregressive modeling approach to the analysis of electronic diary data. BMC Medical Research Methodology, 10(1), 28. doi: 10.1186/1471-2288-10-28."
                       ),
    "GGMncv" = c(
      "Williams, D. (2020). Beyond Lasso: A Survey of Nonconvex Regularization in Gaussian Graphical Models. PsyArXiv pre-print. https://doi.org/10.31234/osf.io/ad57p"
    )
  )
  
  citation <- c(citation,"Epskamp, S., Borsboom, D., & Fried, E. I. (2016). Estimating psychological networks and their accuracy: a tutorial paper. arXiv preprint, arXiv:1604.08462.")
  
  citation
}

print.bootnet <- function(x, ...){
  if (is.list(x$sample$graph)){
    cat("=== bootnet Results (multiple graphs)===")
    print(x$sample)
    cat("\n")
    cat("\nNumber of bootstrapped networks:",length(x[['boots']]),
        paste0("\nResults of original samples stored in ",name,"$sample"),
        paste0("\nTable of all statistics from original samples stored in ",name,"$sampleTable"),
        paste0("\nResults of bootstraps stored in ",name,"$boots"),
        paste0("\nTable of all statistics from bootstraps stored in ",name,"$bootTable"),
        "\n",
        paste0("\nUse plot(",name,"$sample, graph = '...') to plot estimated network of original sample"),
        paste0("\nUse summary(",name,", graph = '...') to inspect summarized statistics (see ?summary.bootnet for details)"),
        paste0("\nUse plot(",name,", graph = '...') to plot summarized statistics (see ?plot.bootnet for details)"),
        "\n\nRelevant references:\n\n",paste0("\t",getRefs(x$sample$default),collapse="\n")
    )
  } else {
    
    directed <- x$sample$directed
    
    if (directed){
      ind <- matrix(TRUE,ncol(x$sample$graph),ncol(x$sample$graph))
    } else {
      ind <- upper.tri(x$sample$graph,diag=FALSE)
    }
    
    
    name <- deparse(substitute(x))[[1]]
    if (nchar(name) > 10) name <- "object"
    cat("=== bootnet Results ===")
    cat("\nNumber of nodes:",nrow(x$sample[['graph']]),
        "\nNumber of non-zero edges in sample:",sum(x$sample[['graph']][ind]!=0),"/",sum(ind),
        "\nMean weight of sample:",mean(x$sample[['graph']][ind]) ,
        "\nNumber of bootstrapped networks:",length(x[['boots']]),
        paste0("\nResults of original sample stored in ",name,"$sample"),
        paste0("\nTable of all statistics from original sample stored in ",name,"$sampleTable"),
        paste0("\nResults of bootstraps stored in ",name,"$boots"),
        paste0("\nTable of all statistics from bootstraps stored in ",name,"$bootTable"),
        "\n",
        paste0("\nUse plot(",name,"$sample) to plot estimated network of original sample"),
        paste0("\nUse summary(",name,") to inspect summarized statistics (see ?summary.bootnet for details)"),
        paste0("\nUse plot(",name,") to plot summarized statistics (see ?plot.bootnet for details)"),
        "\n\nRelevant references:\n\n",paste0("\t",getRefs(x$sample$default),collapse="\n")
    )
  }

}

print.bootnetResult <- function(x, ...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  
  # Trick for printing multiple networks:
  if (is.list(x$graph)){
    cat(paste0("\n=== Estimated networks ==="))
    cat(paste0("\nDefault set used: ",x$default),     
        "\n",
        paste0("\nUse bootnet(",name,") to bootstrap edge weights and centrality indices"),
        "\n\nRelevant references:\n\n",paste0("\t",getRefs(x$default),collapse="\n"))
      
    for (i in 1:length(x$graph)){
      if (x$directed[[i]]){
        ind <- matrix(TRUE,ncol(x$graph[[i]]),ncol(x$graph[[i]]))
      } else {
        ind <- upper.tri(x$graph[[i]],diag=FALSE)
      }
      
      cat(paste0("\n\n=== ",names(x$graph)[[i]]," ==="))
      cat("\nNumber of nodes:",nrow(x[['graph']][[i]]),
          "\nNumber of non-zero edges:",sum(x[['graph']][[i]][ind]!=0),"/",sum(ind),
          "\nMean weight:",mean(x[['graph']][[i]][ind]) ,
          paste0("\nNetwork stored in ",name,"$graph$",names(x$graph)[[i]]),
          paste0("\nUse plot(",name,", graph = '",names(x$graph)[[i]],"') to plot estimated network")
      )
    }
  } else {
    
    
    directed <- x$directed
    
    if (directed){
      ind <- matrix(TRUE,ncol(x$graph),ncol(x$graph))
    } else {
      ind <- upper.tri(x$graph,diag=FALSE)
    }
    
    name <- deparse(substitute(x))[[1]]
    if (nchar(name) > 10) name <- "object"
    cat(paste0("\n=== Estimated network ==="))
    if (isTRUE(x$thresholded)){
      cat("\nNote: network has been thresholded using 'bootThreshold'")
    }
    cat("\nNumber of nodes:",nrow(x[['graph']]),
        "\nNumber of non-zero edges:",sum(x[['graph']][ind]!=0),"/",sum(ind),
        "\nMean weight:",mean(x[['graph']][ind]) ,
        paste0("\nNetwork stored in ",name,"$graph"),
        "\n",
        paste0("\nDefault set used: ",x$default),     
        "\n",
        paste0("\nUse plot(",name,") to plot estimated network"),
        paste0("\nUse bootnet(",name,") to bootstrap edge weights and centrality indices"),
        "\n\nRelevant references:\n\n",paste0("\t",getRefs(x$default),collapse="\n")
    )
  }
  
}