\name{plot.bootnetResult}
\alias{plot.bootnetResult}
\title{
Plot method for bootnetResult objects
}
\description{
Plots the graph using the qgraph package and the \code{\link[qgraph]{qgraph}} function. Defined as \code{qgraph::qgraph(x[['graph']],labels=x[['labels']],...)}
}
\usage{
\method{plot}{bootnetResult}(x, ...)
}
\arguments{
  \item{x}{
A \code{bootnetResult} object
}
  \item{\dots}{
Arguments sent to \code{\link[qgraph]{qgraph}}.
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}