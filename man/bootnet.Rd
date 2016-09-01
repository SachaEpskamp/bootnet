\name{bootnet}
\alias{bootnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrapped network estimation
}
\description{
This function can be used to bootstrap network estimation methods so that the spread of parameter and centrality estimates can be assessed. Most important methods are \code{type = 'nonparametric'} for the non-parametric bootstrap and \code{type = 'case'} for the case-dropping bootstrap. See also Epskamp, Borsboom and Fried (2016) for more details.
}
\usage{
bootnet(data, nBoots = 1000, 
    default = c("none", "EBICglasso", "pcor","IsingFit","IsingLL"), 
    type = c("nonparametric","parametric","node","person","jackknife"),
    model = c("detect","GGM","Ising"), prepFun, prepArgs = list(), estFun, estArgs, 
    graphFun, graphArgs, intFun, intArgs, verbose = TRUE, labels, alpha = 1, 
    subNodes = 2:(ncol(data)-1), subPersons = round(seq(0.25,0.95,length=10)*nrow(data)), 
    computeCentrality = TRUE, propBoot = 1, 
    replacement = TRUE, edgeResample = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame or matrix containing the raw data. Must be numeric, integer or ordered factors.
}
  \item{nBoots}{
Number of bootstraps
}
  \item{default}{
A string indicating the method to use. See documentation at \code{\link[bootnet]{estimateNetwork}}.
}
\item{type}{
 The kind of bootstrap method to use.
}
\item{model}{
 The modeling framework to use. Automatically detects if data is binary or not.
}
  \item{prepFun}{
See documentation at \code{\link[bootnet]{estimateNetwork}}.
}
  \item{prepArgs}{
See documentation at \code{\link[bootnet]{estimateNetwork}}.
}
  \item{estFun}{
See documentation at \code{\link[bootnet]{estimateNetwork}}.
}
  \item{estArgs}{
See documentation at \code{\link[bootnet]{estimateNetwork}}.
}
  \item{graphFun}{
See documentation at \code{\link[bootnet]{estimateNetwork}}.
}
  \item{graphArgs}{
See documentation at \code{\link[bootnet]{estimateNetwork}}.
}
  \item{intFun}{
See documentation at \code{\link[bootnet]{estimateNetwork}}.
}
  \item{intArgs}{
See documentation at \code{\link[bootnet]{estimateNetwork}}.
}
  \item{verbose}{
Logical. Should progress of the function be printed to the console?
}
  \item{labels}{
A character vector containing the node labels. If omitted the column names of the data are used.
}
  \item{alpha}{
  The centrality tuning parameter as used in \code{\link[qgraph]{centrality}}.
}
  \item{subNodes}{
    Range of nodes to sample in node-drop bootstrap
  }
  \item{subPersons}{
  Range of persons to sample in person-drop bootstrap
  }
  \item{computeCentrality}{
    Logical, should centrality be computed?  
  }
  \item{propBoot}{
    Proportion of persons to sample in bootstraps. Set to lower than 1 for m out of n bootstrap
  }
    \item{replacement}{
      Logical, should replacement be used in bootstrap sampling?
  }
  \item{edgeResample}{
      Logical, if \code{TRUE} only the order of edges based on bootstraps is retained and parameter values of the original network are mapped to this new order.
  }
}

\value{
A \code{bootnet} object with the following elements:
\item{sampleTable}{ A data frame containing all estimated values on the real sample. }
\item{bootTable}{ A data frame containing all estimated values on all bootstrapped samples. }
\item{sample}{ A \code{bootnetResult} object with plot and print method containing the estimated network of the real sample.  }
\item{boots}{A list of \code{bootnetResult} objects containing the raw bootstrap results.}
}

\references{
Epskamp, S., Borsboom, D., & Fried, E. I. (2016). Estimating psychological networks and their accuracy: a tutorial paper. arXiv preprint, arXiv:1604.08462.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\seealso{
\code{\link{estimateNetwork}}, \code{\link{differenceTest}}, \code{\link{corStability}}, \code{\link{plot.bootnet}}, \code{\link{summary.bootnet}}
}
\examples{
# BFI Extraversion data from psych package:
library("psych")
data(bfi)
bfiSub <- bfi[,1:25]

# Estimate network:
Network <- estimateNetwork(bfiSub, default = "EBICglasso")

# Centrality indices:
library("qgraph")
centralityPlot(Network)

\dontrun{
# Estimated network:
plot(Network, layout = 'spring') 

### Non-parametric bootstrap ###
# Bootstrap 1000 values, using 8 cores:
Results1 <- bootnet(Network, nBoots = 1000, nCores = 8)

# Plot bootstrapped edge CIs:
plot(Results1, labels = FALSE, order = "sample")

# Plot significant differences (alpha = 0.05) of edges:
plot(Results1, "edge", plot = "difference",onlyNonZero = TRUE,
     order = "sample") 

# Plot significant differences (alpha = 0.05) of node strength:
plot(Results1, "strength", plot = "difference")

# Test for difference in strength between node "A1" and "C2":
differenceTest(Results,  "A1", "C2", "strength")

### Case-drop bootstrap ###
# Bootstrap 1000 values, using 8 cores:
Results2 <- bootnet(Network, nBoots = 1000, nCores = 8, 
              type = "case")

# Plot centrality stability:
plot(Results2)

# Compute CS-coefficients:
corStability(Results2)
}
}
