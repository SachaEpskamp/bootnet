\name{plot.bootnetResult}
\alias{plot.bootnetResult}
\title{
Plot method for bootnetResult objects
}
\description{
Plots the graph using the qgraph package and the \code{\link[qgraph]{qgraph}} function. Defined as \code{qgraph::qgraph(x[['graph']],labels=x[['labels']],...)}
}
\usage{
\method{plot}{bootnetResult}(x, graph, weighted, signed, directed, labels, 
      layout = "spring", parallelEdge = TRUE, cut = 0, theme = "colorblind", ...)
}
\arguments{
  \item{x}{
A \code{bootnetResult} object
}
\item{graph}{Numeric or string indicating which graph to plot. Only needed when multiple graphs are estimated. For example, when using \code{default = "graphicalVAR"}, \code{graph = "temporal"} plots the temporal network and \code{graph = "contemporaneous"} plots the contemporaneous network.}
\item{weighted}{Logical, should the analyzed network be weighted?}
\item{signed}{Logical, should the analyzed network be signed?}
\item{directed}{Logical, is the analyzed network directed? Usually does not have to be set and is detected automatically.}
\item{labels}{Labels of the nodes. Defaults to the column names of the data if missing.}
\item{layout}{Placement of the nodes. See  \code{\link[qgraph]{qgraph}}. Always defaults to \code{"spring"}.}
\item{parallelEdge}{Should edges in directed networks be plotted parallel? See  \code{\link[qgraph]{qgraph}}. Defaults to \code{TRUE} instead of \code{FALSE} (as in qgraph).}
\item{cut}{Should scaling in width and saturation of edges be split? See  \code{\link[qgraph]{qgraph}}. Defaults to \code{0} to disable cut (qgraph chooses a cutoff with at least 20 nodes).}
\item{theme}{Theme of the edge and node colors. See  \code{\link[qgraph]{qgraph}}. Defaults to \code{"colorblind"} rather than the default used in qgraph (\code{"classic"}).}
  \item{\dots}{
Arguments sent to \code{\link[qgraph]{qgraph}}.
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}