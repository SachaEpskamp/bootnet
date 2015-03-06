\name{summary.bootnet}
\alias{summary.bootnet}
\title{
Summarize bootnet results
}
\description{
Creates a data frame (wrapped as \code{\link[dplyr]{tbl_df}}) containing summarized results of the bootstraps.
}
\usage{
\method{summary}{bootnet}(object, types = c("edge", "intercept", "strength", 
                            "closeness", "betweenness", "distance"), ...)
}
\arguments{
  \item{object}{
A \code{bootnet} object
}
  \item{types}{
The types of statistics to include in the summary table
}
\item{\dots}{Not used.}
}
\value{
A \code{tbl_df} (data frame) containing summarized statistics.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}