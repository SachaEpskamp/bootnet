\name{plot.netSimulator}
\alias{plot.netSimulator}
\alias{print.netSimulator}
\alias{summary.netSimulator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
netSimulator S3 methods
}
\description{
Plot, print and summary methods for \code{\link{netSimulator}} output.
}
\usage{
\method{plot}{netSimulator}(x, xvar = "factor(nCases)",
                              yvar = c("sensitivity", "specificity", "correlation"),
                              xfacet = "measure", yfacet = ".", color = NULL,
                             ylim = c(0,1), print = TRUE,  xlab = "Number of cases", 
                              ylab, ...)
\method{print}{netSimulator}(x, digits = 2, ...)
\method{summary}{netSimulator}(object, digits = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Output of \code{\link{netSimulator}}.
}
  \item{object}{
Output of \code{\link{netSimulator}}.
}
  \item{xvar}{
String indicating the variable to be used on the x-axis.
}
  \item{yvar}{
String vector indicating the variable(s) to be used on the y-axis.
}
  \item{xfacet}{
String indicating the variable to be used on the horizontal facets (or \code{"."} to omit).
}
  \item{yfacet}{
String indicating the variable to be used on the vertical facets (or \code{"."} to omit).
}
  \item{color}{
String indicating the variable to be used in coloring boxplots.
}
  \item{ylim}{
Y-axis limits.
}
  \item{print}{
Logical, should the plot be printed? This helps printing the plots to PDF files.
}
\item{digits}{
Number of digits to be used in print and summary method.
}
  \item{xlab}{
X-axis label
}
  \item{ylab}{
Y-axis label. Defaults to combining the values in \code{yvar}. Is hidden when \code{xfacet = "measure"}, as then it is clear what the y-axis represent from the facet labels.
}
  \item{\dots}{
Not used.
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}