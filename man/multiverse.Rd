\name{multiverse}
\alias{multiverse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiverse plot of bootnet results
}
\description{
This function makes a 'multiverse' plot of bootstrap results. Every row indicates an edge and every column a bootstrap; colors are in line of the edge strength as drawn with \code{plot.bootnetResult}.
}
\usage{
multiverse(x, labels = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Results from \code{\link{bootnet}}}
  \item{labels}{
Logical, should labels be printed next to the plot?
}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

