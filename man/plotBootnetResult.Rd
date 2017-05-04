\name{plot.bootnetResult}
\alias{plot.bootnetResult}
\title{
Plot method for bootnetResult objects
}
\description{
Plots the graph using the qgraph package and the \code{\link[qgraph]{qgraph}} function. Defined as \code{qgraph::qgraph(x[['graph']],labels=x[['labels']],...)}
}
\usage{
\method{plot}{bootnetResult}(x, weighted, signed, directed, ...)
}
\arguments{
  \item{x}{
A \code{bootnetResult} object
}
\item{weighted}{Logical, should the analyzed network be weighted?}
\item{signed}{Logical, should the analyzed network be signed?}
\item{directed}{Logical, is the analyzed network directed? Usually does not have to be set and is detected automatically.}
  \item{\dots}{
Arguments sent to \code{\link[qgraph]{qgraph}}.
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}