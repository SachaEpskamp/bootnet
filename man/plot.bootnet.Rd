\name{plot.bootnet}
\alias{plot.bootnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots bootnet results
}
\description{
This function can be used to plot bootnet results by plotting all bootstrapped statistics as line or by plotting confidence intervals (not yet implemented).
}
\usage{
\method{plot}{bootnet}(x, types = c("strength", "closeness", "betweenness"), 
        plot = c("line","interval"), sampleColor = "darkred", samplelwd = 1.1, 
        bootColor = "black", bootAlpha = 0.01, bootlwd = 1, order = c("id","value"),
        decreasing = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{bootnet} object
}
  \item{types}{
The types of statistics to plot. \code{c("intercept", "strength", "closeness", "betweenness") } plots all nodewise statistics whereas \code{c("edge","distance") } plots all pairwise statistics.
}
\item{plot}{
Character string indicating what to plot. \code{"line"} will produce a line graph and \code{"interval"} a 95\% confidence interval graph
}
  \item{sampleColor}{
Color of the original sample line
}
  \item{samplelwd}{
Line width of the original sample line
}
  \item{bootColor}{
Color of the bootstrap lines
}
  \item{bootAlpha}{
Alpha of the bootstrap lines
}
  \item{bootlwd}{
Line width of the bootstrap lines
}
  \item{order}{
    String indicating how to order nodes. \code{"id"} will order nodes based on their name and \code{"value"} will order nodes based on the first entry in \code{type}.
  }
  \item{decreasing}{
    Logical indicating if the ordering is decreasing or increasing.
  }
  \item{\dots}{Not used.}
}

\value{
A \code{ggplot2} object.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
