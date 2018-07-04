\name{genGGM}
\alias{genGGM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a GGM small-world network.
}
\description{
Simulates a GGM as described by Yin and Li (2011), using the Watts and Strogatz (1998) algorithm for generating the graph structure (see \code{\link[igraph]{watts.strogatz.game}}).
}
\usage{
genGGM(Nvar, p = 0, nei = 1, parRange = c(0.5,1), constant = 1.5, propPositive = 0.5,
  clusters = NULL, graph = c("smallworld","random", "scalefree", "hub", "cluster"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Nvar}{
Number of nodes
}
  \item{p}{
Rewiring probability if graph = "smallworld" or "cluster", or connection probability if graph = "random". If cluster, can add multiple p's for each cluster, e.g., "c(.1, .5)" 
}
  \item{nei}{
Neighborhood (see \code{\link[igraph]{watts.strogatz.game}}).
}
  \item{parRange}{
Range of partial correlation coefficients to be originally sampled.
}
  \item{constant}{
A constant as described by  Yin and Li (2011).
}
  \item{propPositive}{
Proportion of edges to be set positive.
}
\item{clusters}{
Number of clusters if graph = "cluster"
}
\item{graph}{Type of graph to simulate}
}
\references{
Yin, J., and Li, H. (2011). A sparse conditional gaussian graphical model for analysis of genetical genomics data. The annals of applied statistics, 5(4), 2630.

Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of `small-world' networks. nature, 393(6684), 440-442.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

