# bootnetResult methods:
print.bootnetResult <- function(x, ...){
  print(x[['graph']])
}

summary.bootnetResult <- function(object, ...){
  cat("\nNumber of nodes:",nrow(object[['graph']]),
      "\nNumber of non-zero edges:",sum(object[['graph']][upper.tri(object[['graph']],diag=FALSE)]==0) ,
      "\nSparsity:",mean(object[['graph']][upper.tri(object[['graph']],diag=FALSE)]) ,
      "\nNumber of intercepts:",NROW(object[['intercepts']])
      )
}

plot.bootnetResult <- function(x,...){
  qgraph::qgraph(x[['graph']],labels=x[['labels']],...)
}
