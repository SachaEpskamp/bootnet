\name{IsingGenerator}
\alias{IsingGenerator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a function that simulates data from the Ising model
}
\description{
Uses \code{\link[IsingSampler]{IsingSampler}} to generate the data.
}
\usage{
IsingGenerator(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Arguments passed to \code{\link[IsingSampler]{IsingSampler}}
}
}
\value{
A function with as first argument the sample size and as second argument a named list, with element \code{graph} encoding a weights matrix and element \code{intercepts} encoding thresholds.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}