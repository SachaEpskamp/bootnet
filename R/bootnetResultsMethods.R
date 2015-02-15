# bootnetResult methods:
print.bootnetResult <- function(x){
  print(x[['graph']])
}

summary.bootnetResult <- function(x){
  cat("\nNumber of nodes:",nrow(x[['graph']]),
      "\nNumber of non-zero edges:",sum(x[['graph']][upper.tri(x[['graph']],diag=FALSE)]==0) ,
      "\nSparsity:",mean(x[['graph']][upper.tri(x[['graph']],diag=FALSE)]) ,
      "\nNumber of intercepts:",NROW(x[['intercepts']])
      )
}

plot.bootnetResult <- function(x,...){
  qgraph::qgraph(x[['graph']],labels=x[['labels']],...)
}
