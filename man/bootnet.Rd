\name{bootnet}
\alias{bootnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrapped network estimation
}
\description{
This function can be used to bootstrap network estimation methods so that the spread of parameter and centrality estimates can be assessed. 
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
A string indicating the method to use. Specifying a \code{default} sets default values to \code{prepFun}, \code{prepArgs}, \code{estFun}, \code{estArgs}, \code{graphFun}, \code{graphArgs}, \code{intFun} and \code{intArgs}. Setting a \code{default} can be omitted but that does require specifying all above mentioned arguments. Current options are:
\describe{
\item{\code{"EBICglasso"}}{Gaussian Markov random field estimation using graphical LASSO and extended Bayesian information criterium to select optimal regularization parameter. Using \code{\link[qgraph]{EBICglasso}} from the qgraph package.}
\item{\code{"IsingFit"}}{Ising model estimation using LASSO regularized nodewise logistic regression and extended Bayesian information criterium to select optimal regularization parameter. Using \code{\link[IsingFit]{IsingFit}} from the IsingFit package.}
\item{\code{"pcor"}}{Partial correlation network (non-regularized Gaussian Markov random field), using \code{\link[corpcor]{cor2pcor}} from the corpcor package.}
}
See details section for a more detailed description.
}
\item{type}{
 The kind of bootstrap method to use.
}
\item{model}{
 The modeling framework to use. Automatically detects if data is binary or not.
}
  \item{prepFun}{
A function that takes as input the raw data and returns whatever the estimation function needs (as first argument). Typically this function is used to correlate or binarize the data. Defaults to \code{\link{identity}} if omitted and \code{default} is not set.
}
  \item{prepArgs}{
A list with arguments for \code{prepFun}. Defaults to \code{list()} if omitted and \code{default} is not set.
}
  \item{estFun}{
A function that takes as input the prepared raw data and estimates a network.
}
  \item{estArgs}{
A list with arguments for \code{estArgs}
}
  \item{graphFun}{
A function that takes the result of \code{estFun} and extracts the estimated weights matrix. Defaults to \code{\link{identity}} if omitted and \code{default} is not set.
}
  \item{graphArgs}{
A list with arguments for \code{graphFun}. Defaults to \code{list()} if omitted and \code{default} is not set.
}
  \item{intFun}{
A function that takes the result of \code{estFun} and extracts the estimated intercepts. Defaults to \code{\link{null}} if omitted and \code{default} is not set.
}
  \item{intArgs}{
A list with arguments for \code{intFun}. Defaults to \code{list()} if omitted and \code{default} is not set.
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

\details{
The following defaults can be used:
\describe{
\item{\code{default = "EBICglasso"}}{
\describe{
\item{prepFun}{\code{\link[qgraph]{cor_auto}} from the qgraph package}
\item{prepArgs}{\code{list(verbose=FALSE)}}
\item{estFun}{\code{\link[qgraph]{EBICglasso}} from the qgraph package}
\item{estArgs}{\code{list(n = nrow(data), returnAllResults = TRUE)}}
\item{graphFun}{\code{function(x)x[['optnet']]}}
\item{graphArgs}{\code{list()}}
\item{intFun}{\code{\link[bootnet]{null}}}
\item{intArgs}{\code{list()}}
}
}
\item{\code{default = "IsingFit"}}{
\describe{
\item{prepFun}{\code{\link[bootnet]{binarize}}}
\item{prepArgs}{\code{list()}}
\item{estFun}{\code{\link[IsingFit]{IsingFit}} from the IsingFit package}
\item{estArgs}{\code{list(plot = FALSE, progress = FALSE)}}
\item{graphFun}{\code{function(x)x[['weiadj']]}}
\item{graphArgs}{\code{list()}}
\item{intFun}{\code{function(x)x[['thresholds']]}}
\item{intArgs}{\code{list()}}
}
}

\item{\code{default = "pcor"}}{
\describe{
\item{prepFun}{\code{\link[qgraph]{cor_auto}} from the qgraph package}
\item{prepArgs}{\code{list()}}
\item{estFun}{\code{\link[corpcor]{cor2pcor}} from the corpcor package}
\item{estArgs}{\code{list()}}
\item{graphFun}{\code{\link{identity}}}
\item{graphArgs}{\code{list()}}
\item{intFun}{null}
\item{intArgs}{\code{list()}}
}
}
}
}
\value{
A \code{bootnet} object with the following elements:
\item{sampleTable}{ A data frame containing all estimated values on the real sample. }
\item{bootTable}{ A data frame containing all estimated values on all bootstrapped samples. }
\item{sample}{ A \code{bootnetResult} object with plot and print method containing the estimated network of the real sample.  }
\item{boots}{A list of \code{bootnetResult} objects containing the raw bootstrap results.}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com
}
\seealso{
\code{\link{plot.bootnet}}, \code{\link{summary.bootnet}}
}
\examples{
# BFI Extraversion data from psych package:
library("psych")
data(bfi)
bfiSub <- bfi[,11:15]

# 10 bootstraps (using Pearson correlations to speed up):
Results <- bootnet(bfiSub, nBoot = 5, default = "EBICglasso", prepFun = cor, 
                   prepArgs = list(use = "pairwise.complete.obs"))


# Results:
print(Results)

# Estimated network:
plot(Results$sample, layout = 'spring') 

# Estimated values in sample:
Results$sampleTable

# Estimated samples in bootstraps:
Results$bootTable

# Summary table:
summary(Results)

# Summary of strength only:
summary(Results, "strength")

# Plots of edges:
plot(Results, plot = "area", order = "sample")
plot(Results, plot = "interval", order = "sample")

#### Centrality ####

# Node-drop:
# 100 bootstraps (using Pearson correlations to speed up):
NodeDrop <- bootnet(bfiSub, nBoot = 25, default = "EBICglasso", prepFun = cor, 
                   prepArgs = list(use = "pairwise.complete.obs"), type = "node")

# Person-drop:
PersonDrop <- bootnet(bfiSub, nBoot = 25, default = "EBICglasso", prepFun = cor, 
                    prepArgs = list(use = "pairwise.complete.obs"), type = "node")

\dontrun{
  # 1000 bootstraps (using Pearson correlations to speed up):
  NodeDrop <- bootnet(bfiSub, nBoot = 1000, default = "EBICglasso", prepFun = cor, 
                      prepArgs = list(use = "pairwise.complete.obs"), type = "node")
  
  # Person-drop:
  PersonDrop <- bootnet(bfiSub, nBoot = 1000, default = "EBICglasso", prepFun = cor, 
                        prepArgs = list(use = "pairwise.complete.obs"), type = "node")
}

# Plot:
plot(NodeDrop) # Correlation with original
plot(NodeDrop, perNode = TRUE) # Value per sampling level

plot(PersonDrop) # Correlation with original
plot(PersonDrop, perNode = TRUE) # Value per sampling level
}
