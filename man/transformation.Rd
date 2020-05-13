\name{transformation}
\alias{quantile_transformation}
\alias{rank_transformation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data transformation functions
}
\description{
Functions to transform data
}
\usage{
quantile_transformation(x)
rank_transformation(x, ties.method = c("average", "first", 
    "last", "random", "max", "min"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A dataset
}
\item{ties.method}{See \code{rank}.}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
