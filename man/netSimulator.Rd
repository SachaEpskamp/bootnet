\name{netSimulator}
\alias{netSimulator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Network Estimation Performance
}
\description{
This function can be used to run a simulation study on the performance of network estimation given various sample sizes. The purpose of this function is to provide a way to assess the required sample size given a network structure. By default, the function uses \code{\link{simGraph}} to simulate a chain graph or small-world network. See details for more information.
}
\usage{
netSimulator(input = simGraph(Nvar = 10), dataGenerator = ggmGenerator(ordinal = FALSE), nCases = c(50, 100, 250, 500, 1000, 2500), nReps = 100, nCores = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
Either a weights matrix, a list containing elements \code{graph} (encoding the weights matrix) and \code{intercepts} (encoding the intercepts), or a function generating such objects.
}
  \item{dataGenerator}{
A function that generates data. The first argument must be the sample size, the second argument must be the output of \code{input}.
}
  \item{nCases}{
The sample sizes to test for.
}
  \item{nReps}{
Number of repetitions per sampling level.
}
  \item{nCores}{
Number of cores to use. Set to more than 1 to use parallel computing.
}
  \item{\dots}{
Arguments used by \code{\link{estimateNetwork}} to estimate the network structure.
}
}
\details{
The \code{dataGenerator} argument can be any function that generates data. Currently, only \code{\link{ggmGenerator}} is implemented in bootnet itself, which generates data given a Gaussian graphical model.
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\examples{
# 5-node GGM chain graph:
Res <- netSimulator(simGraph(5), nReps = 10)
Res
plot(Res)

\dontrun{
# BFI example:
# Load data:
library("psych")
data(bfi)
bfiData <- bfi[,1:25]

# Estimate a network structure, with parameters refitted without LASSO regularization:
library("qgraph")
Network <- EBICglasso(cor_auto(bfiData), nrow(bfiData), refit = TRUE)

# Treat the estimated network as the true network:
trueNetwork <- Network

# Run simulation study
# 
# - input: input to the generator, can be a weights matrix or generating function
# - dataGenerator: A function that generates data given a network, ggmGenerator creates such a function.
# - nCores: Number of computer cores
# - nReps: Number of repetitions per condition
# - nCases: Different conditions for number of cases to test.
library("bootnet")
simResults <- netSimulator(input = trueNetwork,
                           dataGenerator = ggmGenerator(ordinal=TRUE,
                                                        nLevels = 6),
                           nCases = c(100,250,500,1000),
                           nReps = 100,
                           nCores = 6,
                           default = "EBICglasso"
                           )

# Average values:
simResults

# Plot boxplots:
plot(simResults)
}
}
