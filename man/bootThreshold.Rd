\name{bootThreshold}
\alias{bootThreshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Threshold network based on bootstrapped intervals
}
\description{
This function takes the output of \code{\link{bootnet}} and returns a network as if it had been estimated using \code{\link{estimateNetwork}}, but with edges removed (set to zero) based on some significance level.
}
\usage{
bootThreshold(bootobject, alpha = 0.05, verbose = TRUE, thresholdIntercepts = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bootobject}{
Nonparametric bootstrap results from \code{\link{bootnet}}
}
  \item{alpha}{
Significance level
}
  \item{verbose}{
Logical, should progress be reported to the console?
}
  \item{thresholdIntercepts}{
Logical, should intercepts also be thresholded?
}
}
\value{
A \code{bootnetResult} object with the following elements:
\item{graph}{The weights matrix of the network}
\item{intercepts}{The intercepts}
\item{results}{The results of the estimation procedure}
\item{labels}{A vector with node labels}
\item{nNodes}{Number of nodes in the network}
\item{nPerson}{Number of persons in the network}
\item{input}{Input used, including the result of the default set used}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\seealso{
\code{\link{bootnet}}, \code{\link{estimateNetwork}}
}
\examples{
\dontrun{
# BFI Extraversion data from psych package:
library("psych")
data(bfi)
bfiSub <- bfi[,1:25]

# Estimate unregularized network:
Network <- estimateNetwork(bfiSub, default = "pcor", corMethod = "cor")

# Bootstrap 1000 values, using 8 cores:
boots <- bootnet(Network, nBoots = 1000, nCores = 8)

# Threshold network:
Network_thresholded <- bootThreshold(boots)

# Plot:
plot(Network_thresholded)
}
}

