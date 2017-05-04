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
bootnet(data, nBoots = 1000, default = c("none", "EBICglasso", "pcor","IsingFit",
                  "IsingSampler", "huge","adalasso","mgm"), type =
                 c("nonparametric", "parametric", "node", "person",
                 "jackknife", "case"), nCores = 1, statistics = c("edge"
                 ,"strength","closeness","betweenness"), model = c("detect",
                 "GGM", "Ising"), fun, prepFun, prepArgs, estFun, estArgs,
                 graphFun, graphArgs, intFun, intArgs, verbose = TRUE,
                 construct = c("default", "function", "arguments"), labels, 
                 alpha = 1, caseMin = 0.05, caseMax = 0.75, caseN = 10,
                 subNodes = 2:(ncol(data) - 1), subCases = round((1 -
                 seq(caseMin, caseMax, length = caseN)) * nrow(data)),
                 computeCentrality = TRUE, propBoot = 1, replacement =
                 TRUE, graph, sampleSize, intercepts, weighted, signed,
                 directed, ...)
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
\item{nCores}{
Number of cores to use in computing results. Set to 1 to not use parallel computing.
}
\item{statistics}{
Vector indicating which statistics to store. Can contain \code{"edge"}, \code{"strength"}, \code{"closeness"}, \code{"betweenness"}, \code{"length"} and \code{"distance"}. By default, length and distance are not stored.
}
\item{model}{
 The modeling framework to use. Automatically detects if data is binary or not.
}
\item{fun}{
A custom estimation function, when no default set is used. This must be a function that takes the data as input (first argument) and returns either a weights matrix or a list containing the elements \code{"graph"} for the weights matrix, \code{"intercepts"} for the intercepts (optional) and \code{"results"} for the full estimation results (optional).
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
  \item{caseMin}{Minimum proportion of cases to drop when \code{type = "case"}.}
    \item{caseMax}{Maximum proportion of cases to drop when \code{type = "case"}.}
      \item{caseN}{Number of sampling levels to test when \code{type = "case"}.}
  \item{subCases}{
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
  \item{graph}{A given network structure to use in parametric bootstrap.}
  \item{sampleSize}{The samplesize to use in parametric bootstrap.}
  \item{intercepts}{Intercepts to use in parametric bootstrap.}
  \item{construct}{How to construct the estimator? Set to \code{"arguments"} to use old functionality (not recommended).}
\item{weighted}{Logical, should the analyzed network be weighted?}
\item{signed}{Logical, should the analyzed network be signed?}
\item{directed}{Logical, is the analyzed network directed? Usually does not have to be set and is detected automatically.}
  \item{\dots}{
  Additional arguments used in the estimator function.
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
