\name{summary.bootnet}
\alias{summary.bootnet}
\title{
Summarize bootnet results
}
\description{
Creates a data frame (wrapped as \code{\link[dplyr]{tbl_df}}) containing summarized results of the bootstraps.
}
\usage{
\method{summary}{bootnet}(object, graph, statistics = c("edge", "intercept",
                 "strength", "closeness", "betweenness", "distance"),
                 perNode = FALSE, rank = FALSE, tol =
                 sqrt(.Machine$double.eps), ...)
}
\arguments{
  \item{object}{
A \code{bootnet} object
}
\item{graph}{Numeric or string indicating which graph to summarize. Only needed when multiple graphs are estimated. For example, when using \code{default = "graphicalVAR"}, \code{graph = "temporal"} plots the temporal network and \code{graph = "contemporaneous"} plots the contemporaneous network.}
  \item{statistics}{
The types of statistics to include in the summary table
}
  \item{perNode}{
    Logical, should centrality estimates per node be plotted instead of correlation with original parameter. Only used in node and person drop bootstrap.
  }
  \item{rank}{
  Logical, should plots show rank of statistics instead of statistics?
}
\item{tol}{Tolerance level for assuming an edge is set to zero.}
\item{\dots}{Not used.}
}
\value{
A \code{tbl_df} (data frame) containing summarized statistics.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}