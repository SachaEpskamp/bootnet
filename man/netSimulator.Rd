\name{netSimulator}
\alias{netSimulator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Network Estimation Performance
}
\description{
This function can be used to run a simulation study on the performance of network estimation by varying sample size or any argument used as input to \code{\link{estimateNetwork}}. The purpose of this function is to provide a way to assess the required sample size given a network structure, as well as to easily perform simulation studies. By default, the function uses \code{\link{genGGM}} to simulate a chain graph or small-world network. See details for more information.
}
\usage{
netSimulator(
      input = genGGM(Nvar = 10), 
      nCases = c(50, 100, 250, 500, 1000, 2500), 
      nReps = 100,
      nCores = 1, 
      default,
      dataGenerator, 
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
Either a weights matrix, a list containing elements \code{graph} (encoding the weights matrix) and \code{intercepts} (encoding the intercepts), or a function generating such objects. By default, \code{\link{genGGM}} is used to generate a Gaussian graphical model. However, it is reccomended to replace this with a prior expected graph structure.
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
  \item{default}{
Default set used (see \code{\link{estimateNetwork}}). In most cases, this will set \code{dataGenerator} to the relevant generator.
}
  \item{dataGenerator}{
A function that generates data. The first argument must be the sample size, the second argument must be the output of \code{input}. Can often be ignored if \code{default} is set.
}
  \item{\dots}{
Arguments used by \code{\link{estimateNetwork}} to estimate the network structure. Providing a vector for any argument will simulate under each value. This way, any argument in \code{\link{estimateNetwork}} can be used in a simulation study.
}
}
\details{
*any* argument to \code{\link{estimateNetwork}} can be used in a simulation study, with a vector (e.g., \code{rule = c("AND","OR")}) specifying that both conditions are tested. Adding too many conditions can quickly make any simulation study intractible, so only vary some arguments! The \code{dataGenerator} argument can be any function that generates data. Currently, only \code{\link{ggmGenerator}} and \code{\link{IsingGenerator}} are implemented in bootnet itself, which generates data given a Gaussian graphical model.
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\examples{
# 5-node GGM chain graph:
trueNetwork <- genGGM(5)

# Simulate:
Res <- netSimulator(trueNetwork, nReps = 10)

# Results:
Res

# Plot:
plot(Res)

\dontrun{
library("bootnet")

# BFI example:
# Load data:
library("psych")
data(bfi)
bfiData <- bfi[,1:25]

# Estimate a network structure, with parameters refitted without LASSO regularization:
library("qgraph")
Network <- EBICglasso(cor_auto(bfiData), nrow(bfiData), refit = TRUE)

# Simulate 100 repititions in 8 cores under different sampling levels:
Sim1 <- netSimulator(Network,
                     default = "EBICglasso",
                     nCases = c(100,250,500),
                     nReps = 100,
                     nCores = 8)

# Table of results:
Sim1

# Plot results:
plot(Sim1)

# Compare different default set at two sampling levels:
Sim2 <- netSimulator(Network,
                     default = c("EBICglasso","pcor","huge"),
                     nCases = c(100,250,500),
                     nReps = 100,
                     nCores = 8)

# Print results:
Sim2

# Plot results:
plot(Sim2, xfacet = "default", yvar = "correlation")

# Difference using polychoric or pearson correlations in ordinal data:
Sim3 <- netSimulator(Network,
                     dataGenerator = ggmGenerator(ordinal = TRUE, nLevels = 4),
                     default = "EBICglasso",
                     corMethod = c("cor","cor_auto"),
                     nCases = c(100,250, 500),
                     nReps = 100,
                     nCores = 8)

# Print results:
Sim3

# Plot results:
plot(Sim3, color = "corMethod")

# Ising model:
trueNetwork <- read.csv('http://sachaepskamp.com/files/weiadj.csv')[,-1]
trueNetwork <- as.matrix(trueNetwork)
Symptoms <- rownames(trueNetwork) <- colnames(trueNetwork)
Thresholds <- read.csv('http://sachaepskamp.com/files/thr.csv')[,-1]

# Create an input list (intercepts now needed)
input <- list(graph=trueNetwork,intercepts=Thresholds)

# Simulate under different sampling levels:
Sim4 <- netSimulator(
  input = input,
  default = "IsingFit",
  nCases = c(250,500,1000),
  nReps = 100,
  nCores = 8)

# Results:
Sim4

# Plot:
plot(Sim4)

# Compare AND and OR rule:
Sim5 <- netSimulator(
  input = input,
  default = "IsingFit",
  nCases = c(250,500,1000),
  rule = c("AND","OR"),
  nReps = 100,
  nCores = 8)

# Print:
Sim5

# Plot:
plot(Sim5, yfacet = "rule")

}
}
