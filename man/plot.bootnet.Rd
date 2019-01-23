\name{plot.bootnet}
\alias{plot.bootnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots bootnet results
}
\description{
This function can be used to plot bootnet results by plotting all bootstrapped statistics as line or by plotting confidence intervals.
}
\usage{
\method{plot}{bootnet}(x, statistics, plot, graph, CIstyle = c(
                 "quantiles", "SE"), rank = FALSE, sampleColor =
                 "darkred", samplelwd = 1, meanColor = "black", meanlwd
                 = 0.5, bootColor = "black", bootAlpha = 0.01, bootlwd
                 = 0.9, areaAlpha = 0.1, order = c("id", "sample",
                 "mean"), decreasing = TRUE, perNode = FALSE,
                 legendNcol = 2, labels = TRUE, legend = TRUE,
                 subsetRange = c(100, 0), area = !perNode, alpha =
                 0.05, onlyNonZero = FALSE, differenceShowValue,
                 differenceEdgeColor = TRUE, verbose = TRUE, panels =
                 TRUE, split0 = FALSE, prop0 = ifelse(split0, TRUE,
                 FALSE), prop0_cex = 1, prop0_alpha = 0.8,
                 prop0_minAlpha = 0.25, subset, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{bootnet} object
}
\item{statistics}{
The types of statistics to plot.Defaults to \code{"edge"} for regular bootstrap and \code{c("strength", "outStrength", "inStrength")} for node and person drop bootstrap. Use \code{"all"} to obtain all pairwise statistics tested for regular bootstraps and all nodewise statistics tested for person and node drop bootstarps.}
\item{plot}{
Character string indicating what to plot. Can be \code{"area"} to produce a graph with the area indicating the confidence region, or \code{"difference"} producing a plot showing significant differences. Other options are \code{"line"} and \code{"interval"}, which are currently unstable and not recommended to use.
}
\item{graph}{If multipe graphs are estimated, which graph should be plotted? Currently used for \code{default = "graphicalVAR"} to plot a temporal network using \code{graph = "temporal"} or a contemporaneous network using \code{graph = "contemporaneous"}}
\item{CIstyle}{
 Style of CIs to construct. \code{"SE"} shows the sample statistic plus and minus two times the standard deviation of bootstraps, and \code{"quantiles"} the area between the 2.5th and 97.5th quantile. Defaults to  \code{"quantiles"}.
}
\item{rank}{
  Logical, should plots show rank of statistics instead of statistics?
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
\item{areaAlpha}{
Alpha of the area
}
  \item{order}{
    String indicating how to order nodes. \code{"id"} will order nodes based on their name, \code{"mean"} will order nodes based on the average bootstrapped value of the first statistic in \code{statistics}, and \code{"sample"} will order the nodes as done in  \code{"mean"} but orders ties based on their sample value.
  }
  \item{decreasing}{
    Logical indicating if the ordering is decreasing or increasing.
  }
  \item{perNode}{
    Logical, should centrality estimates per node be plotted instead of correlation with original parameter. Only used in node and person drop bootstrap.
  }
  \item{legendNcol}{
  Number of columns in the legend if \code{perNode = TRUE}.
  }
  \item{labels}{
    Logical, should labels be plotted?
  }
  \item{legend}{
    Logical, should the legend be plotted?
  }
  \item{subsetRange}{
  Range in percentages of the x-axis in node and person drop plots.
  }
  \item{area}{
   Logical, should the confidence area be plotted?
  }
    \item{alpha}{
   Signifiance level used in \code{plot = "difference"}.
  }
  \item{onlyNonZero}{Logical used when \code{plot = "difference", statistics = "edge"}, should only edges be included that were nonzero in the estimated network structure?}
  \item{differenceShowValue}{Logical used when \code{plot = "difference"}. Should values be shown in the diagonal of the difference plot?
  }
   \item{differenceEdgeColor}{Logical used when \code{plot = "difference", statistics = "edge"}. Should diagonal blocks be colored according to default edge colors?
  }
  \item{verbose}{Should expected alpha be printed?}
  \item{panels}{Logical, should panel titles be printed?}
  \item{meanColor}{Color of the bootstrap means.}
  \item{meanlwd}{Line width of the bootstrap means}
  \item{split0}{Logical. When set to TRUE, the displayed intervals are based on occasions when the parameter was not estimated to be zero, and an extra box is added indicating the number of times a parameter is estimated to be zero.}
    \item{prop0}{Logical, should boxes indicating the proportion of times parameters were estimated to be zero be added to the plot?}
    \item{prop0_cex}{Size of the boxes indicating number of times a parameter was set to zero.}
    \item{prop0_alpha}{Transparency of the boxes indicating number of times a parameter was set to zero.}
    \item{prop0_minAlpha}{Minimal transparency of the *lines* of plotted intervals as the proportion of times an edge was not included goes to 0.}
    \item{subset}{Vector indicating labels of nodes to include in the plot. This can be used to show, for example, only edges related to one particular node.}
  \item{\dots}{Not used.}
}

\value{
A \code{ggplot2} object.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
